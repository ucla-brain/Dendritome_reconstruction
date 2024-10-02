//
// Created by muyezhu on 4/5/19.
//

#include "mcp3d_image_macros.hpp"
#include "mcp3d_image_maths.hpp"
#include "mcp3d_image_base.hpp"

using namespace std;


int64_t mcp3d::MImageBase::VoxelAddress(int c, int z, int y, int x, int t) const
{
    if (empty())
        return -1;
    return mcp3d::LinearAddressRobust(loaded_xyz_dims(), z, y, x);
}

bool mcp3d::MImageBase::HasEqualData(const MImageBase &other)
{
    if (loaded_n_channels() == 0 || other.loaded_n_channels() == 0)
        return false;
    if (loaded_dims() != other.loaded_dims())
        return false;
    if (loaded_voxel_type() != other.loaded_voxel_type())
        return false;
    for (int i = 0; i < loaded_n_channels(); ++i)
        if (!mcp3d::VolumeEqual(data(i), other.data(i), loaded_xyz_dims(), mcp3d::VoxelSize(loaded_voxel_type())))
            return false;
    return true;
}
