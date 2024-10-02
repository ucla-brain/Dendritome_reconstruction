#include <iostream>
#include <vector>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image/mcp3d_image.hpp"
#include "vaa3d/swc_utility/resample_swc/resample_swc.hpp"
#include "app2_parameters.hpp"
#include "all_path_pruning2.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    mcp3d::App2CommandLineArgs args;
    // if command line arguments invalid, do not execute further
    if (!mcp3d::ParseApp2Args(argc, argv, args))
        return 1;
    // read data
	mcp3d::MImage image(args.image_root_dir());
    try
    {
        image.ReadImageInfo(args.channel_number());
        if (image.image_info().empty())
        {
            MCP3D_MESSAGE("no supported image formats found in " +
                          args.image_root_dir() + ", do nothing.")
            return 0;
        }
        image.SaveImageInfo();
        // use unit strides only
        mcp3d::MImageBlock block(args.image_offsets(), args.image_extents());
        image.SelectView(block, args.channel_number(), args.resolution_level());
        image.ReadData();
    }
    catch(...)
    {
        MCP3D_MESSAGE("error in image io. neuron tracing not performed")
        return 1;
    }
    // if extents are not given from command line, fill in the max extents
    // at selected resolution level from image view
    args.set_image_extents(image.loaded_view().pyr_level_extents(args.channel_number()));
	args.PrintParameters();

    
    bool success = mcp3d::AllPathPruning2(image, args);
    if (success)
        mcp3d::ResampleSwc(args.swc_path(), 10.0);

	return 0;
}
