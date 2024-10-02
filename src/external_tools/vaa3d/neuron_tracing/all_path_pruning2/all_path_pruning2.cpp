//
// Created by muyezhu on 5/24/19.
//
#include <cstdint>
#include <stdexcept>
#include "fastmarching_tree.h"
#include "fastmarching_dt.h"
#include "hierarchy_prune.h"
#include "vaa3d/swc_utility/swc_file.hpp"
#include "vaa3d/swc_utility/markers.h"
#include "app2_image_preprocess.hpp"
#include "all_path_pruning2.hpp"

using namespace std;

bool mcp3d::AllPathPruning2(mcp3d::MImageBase &image, mcp3d::App2CommandLineArgs& args)
{
    if (image.loaded_view().empty())
    {
        MCP3D_MESSAGE("no image data loaded. do nothing")
        return false;
    }

    vector<MyMarker> fast_marching_source;
    try
    {
        fast_marching_source = mcp3d::App2ImageProprocessing(image, args);
    }
    catch (const exception& e)
    {
        cout << e.what() << endl;
        MCP3D_MESSAGE("error in image processing, neuron tracing not performed")
        return false;
    }

    int channel_number = image.loaded_view().view_channels()[0];
    int xdim = image.loaded_view().xdim(channel_number);
    int ydim = image.loaded_view().ydim(channel_number);
    int zdim = image.loaded_view().zdim(channel_number);
    mcp3d::VoxelType vtype = image.loaded_view().voxel_type();
    mcp3d::App2Parameters& app2_parameters = args.app2_parameters();

    float * phi = 0;
    // when no marker files is given, use vaa3d's cell body detection
    if (fast_marching_source.empty())
    {
        cout<<"Start detecting cellbody"<<endl;
        if (vtype == mcp3d::VoxelType::M8U)
            fastmarching_dt_XY<uint8_t>(image.Volume<uint8_t>(channel_number), phi, xdim, ydim, zdim,
                                        app2_parameters.cnn_type(), 1);
        else
            fastmarching_dt_XY<uint16_t>(image.Volume<uint16_t>(channel_number), phi, xdim, ydim, zdim,
                                         app2_parameters.cnn_type(), 1);

        int64_t n_plane = (int64_t)xdim * (int64_t)ydim;
        int64_t n_volume = n_plane * zdim;

        int64_t max_loc = 0;
        double max_val = phi[0];
        for(int64_t i = 0; i < n_volume; i++)
        {
            if(phi[i] > max_val)
            {
                max_val = phi[i];
                max_loc = i;
            }
        }
        MyMarker max_marker(max_loc % xdim, max_loc % n_plane / xdim, max_loc / n_plane);
        fast_marching_source.push_back(max_marker);
    }

    cout<<"======================================="<<endl;
    cout<<"Construct neuron tree"<<endl;
    vector<MyMarker *> outtree;
    if (fast_marching_source.empty())
    {
        MCP3D_MESSAGE("need at least one marker")
        return false;
    }
    else if(fast_marching_source.size() == 1)
    {
        MCP3D_MESSAGE("only one input marker")
        if(app2_parameters.gsdt())
        {
            if (phi == 0)
            {
                cout<<"processing fastmarching distance transformation ..."<<endl;
                if (vtype == mcp3d::VoxelType::M8U)
                    fastmarching_dt(image.Volume<uint8_t>(channel_number), phi, xdim, ydim, zdim,
                                    app2_parameters.cnn_type(), 1);
                else
                    fastmarching_dt(image.Volume<uint16_t>(channel_number), phi, xdim, ydim, zdim,
                                    app2_parameters.cnn_type(), 1);

            }

            cout<<endl<<"constructing fastmarching tree ..."<<endl;
            fastmarching_tree(fast_marching_source[0], phi, outtree, xdim, ydim, zdim,
                              app2_parameters.cnn_type(), 1, app2_parameters.allow_gap());
        }
        else
        {
            if (vtype == mcp3d::VoxelType::M8U)
                fastmarching_tree(fast_marching_source[0], image.Volume<uint8_t>(channel_number), outtree, xdim, ydim, zdim,
                                  app2_parameters.cnn_type(), 1, app2_parameters.allow_gap());

            else
                fastmarching_tree(fast_marching_source[0], image.Volume<uint16_t>(channel_number), outtree, xdim, ydim, zdim,
                                  app2_parameters.cnn_type(), 1, app2_parameters.allow_gap());
        }
    }
    else
    {
        vector<MyMarker> target;
        target.insert(target.end(), fast_marching_source.begin() + 1, fast_marching_source.end());
        if (app2_parameters.gsdt())
        {
            if(phi == 0)
            {
                cout << "processing fastmarching distance transformation ..." << endl;
                if (vtype == mcp3d::VoxelType::M8U)
                    fastmarching_dt(image.Volume<uint8_t>(channel_number), phi, xdim, ydim, zdim,
                                    app2_parameters.cnn_type(), 1);
                else
                    fastmarching_dt(image.Volume<uint16_t>(channel_number), phi, xdim, ydim, zdim,
                                    app2_parameters.cnn_type(), 1);
            }
            cout << endl << "constructing fastmarching tree ..." << endl;
            fastmarching_tree(fast_marching_source[0], target, phi, outtree, xdim, ydim, zdim, app2_parameters.cnn_type());
        }
        else
        {
            if (vtype == mcp3d::VoxelType::M8U)
                fastmarching_tree(fast_marching_source[0], target, image.Volume<uint8_t>(channel_number), outtree, xdim, ydim, zdim, app2_parameters.cnn_type());
            else
                fastmarching_tree(fast_marching_source[0], target, image.Volume<uint16_t>(channel_number), outtree, xdim, ydim, zdim, app2_parameters.cnn_type());
        }
    }

    cout << endl << "Pruning neuron tree" << endl;
    vector<MyMarker*> & inswc = outtree;
    vector<MyMarker*> outswc;
    // always use hierarchical coverage pruning (APP2 paper)
    if (vtype == mcp3d::VoxelType::M8U)
        happ(inswc, outswc, image.Volume<uint8_t>(channel_number), xdim, ydim, zdim,
             1, app2_parameters.length_thresh(), app2_parameters.sr_ratio());
    else
        happ(inswc, outswc, image.Volume<uint16_t>(channel_number), xdim, ydim, zdim,
             1, app2_parameters.length_thresh(), app2_parameters.sr_ratio());

    saveSWC_file(args.swc_path(), outswc, args.MetaString());

    if (phi)
    {
        delete[] phi;
        phi = 0;
    }
    for(size_t i = 0; i < outtree.size(); i++)
        delete outtree[i];
    return true;
}

