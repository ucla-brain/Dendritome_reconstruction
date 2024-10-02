//
// Created by muyezhu on 7/31/19.
//

#include "common/mcp3d_macros.hpp"
#include "image/mcp3d_image_macros.hpp"
#include "app2_image_preprocess.hpp"

vector<MyMarker> mcp3d::App2ImageProprocessing(mcp3d::MImageBase &image,
                                               mcp3d::App2CommandLineArgs &args)
{
    MCP3D_ASSERT(!image.loaded_view().empty())

    if (image.loaded_view().n_channels() > 1)
        MCP3D_MESSAGE("multiple channel data loaded. will only reconstruct from the first of selected channels")
    int channel_number = image.loaded_view().view_channels()[0];

    mcp3d::VoxelType vtype = image.loaded_view().voxel_type();
    if (vtype != mcp3d::VoxelType::M8U && vtype != mcp3d::VoxelType::M16U)
        MCP3D_DOMAIN_ERROR("unsupported voxel type")

    mcp3d::App2Parameters& app2_parameters = args.app2_parameters();
    vector<MyMarker> fast_marching_source;

    int marker_intensity = -1;
    if(!app2_parameters.marker_file_path().empty())
    {
        fast_marching_source = readMarker_file(app2_parameters.marker_file_path());
        // enforcing single input marker (aka soma) for now
        MCP3D_ASSERT(fast_marching_source.size() == 1)
        for (auto& inmarker: fast_marching_source)
        {
            // scale and adjust offset for marker
            // (input markers are global and at pyr level 0)
            vector<int> zyx_scales {image.image_info().pyr_z_ratios()[args.resolution_level()],
                                    image.image_info().pyr_xy_ratios()[args.resolution_level()],
                                    image.image_info().pyr_xy_ratios()[args.resolution_level()]};
            inmarker /= zyx_scales;
            inmarker -= image.loaded_view().pyr_level_offsets(channel_number);
            inmarker.floor_coordinates();

            IMAGE_LOADED_TYPED_CALL(::MarkerIntensity, image, image, inmarker, marker_intensity);
            if (marker_intensity < 0)
                MCP3D_OUT_OF_RANGE("marker coordinates out of image bounds. exiting.")
            cout << "marker at (zyx) " << inmarker.z << ", " << inmarker.y << ", " << inmarker.x << endl;
            cout << "marker voxel intensity = " << marker_intensity << endl;
        }
    }

    // if foreground percentile is given, set corresponding percentage of voxels
    // to foreground
    if (app2_parameters.foreground_percent() > 0)
    {
        // foreground percentile applied to entire volume
        if (!app2_parameters.foreground_percent_by_plane())
        {
            // perform threshold and set app2_parameters.threshold_value
            IMAGE_LOADED_TYPED_CALL(::App2ThresholdVolumeTopPercentile, image,
                                    image, app2_parameters);
        }
            // foreground percentile applied to each xy plane
        else
        {
            // perform threshold and set app2_parameters.plane_threshold_values
            IMAGE_LOADED_TYPED_CALL(::App2ThresholdPlanesTopPercentile, image,
                                    image, app2_parameters);
        }
    }
    // otherwise threshold entire volume with with app2_parameters.threshold_value
    else
    {
        IMAGE_LOADED_TYPED_CALL(::App2ThresholdVolume, image, image, app2_parameters);
    }

    if (!fast_marching_source.empty())
    {
        int thresholded_marker_intensity = -1;
        IMAGE_LOADED_TYPED_CALL(::MarkerIntensity, image, image, fast_marching_source[0], thresholded_marker_intensity);
        if (thresholded_marker_intensity == 0)
        {
            int plane_threshold_value = app2_parameters.foreground_percent_by_plane() ?
                                        app2_parameters.plane_threshold_values()[fast_marching_source[0].z] :
                                        app2_parameters.threshold_value();
                MCP3D_RUNTIME_ERROR("marker voxel intensity less than threshold value = " +
                                    to_string(plane_threshold_value)+ ", exiting")
        }

    }

    if (app2_parameters.foreground_percent() > 0)
    {
        if (!app2_parameters.foreground_percent_by_plane())
            cout << "setting " << app2_parameters.foreground_percent() * 100
                 << "% voxels to foreground, equivalent of threshold value = "
                 << app2_parameters.threshold_value() << endl;
        else
            cout << "setting " << app2_parameters.foreground_percent() * 100
                 << "% voxels to foreground per plane, equivalent of plane threshold values = "
                 << app2_parameters.plane_threshold_values_string() << endl;
    }
    else
        cout << "setting background thresh value = " << app2_parameters.threshold_value() << endl;

    return fast_marching_source;
}

