//
// Created by muyezhu on 7/30/19.
//

#ifndef MCP3D_APP2_IMAGE_PREPROCESS_HPP
#define MCP3D_APP2_IMAGE_PREPROCESS_HPP

#include <limits>
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image/mcp3d_image_maths.hpp"
#include "image/mcp3d_image_base.hpp"
#include "app2_parameters.hpp"
#include "vaa3d/swc_utility/markers.h"

namespace mcp3d
{

std::vector<MyMarker> App2ImageProprocessing(MImageBase &image, App2CommandLineArgs& args);

}

template <typename VType>
void MarkerIntensity(mcp3d::MImageBase &image, MyMarker& marker, int& marker_intensity)
{
    MCP3D_ASSERT(image.n_volumes() == 1)
    int channel_number = image.loaded_view().view_channels()[0];
    if (marker.z < 0 || marker.z >= image.loaded_view().zdim(channel_number) ||
        marker.y < 0 || marker.y >= image.loaded_view().ydim(channel_number) ||
        marker.x < 0 || marker.x >= image.loaded_view().xdim(channel_number))
        return;
    VType intensity = image.operator()<VType>(channel_number, (int)marker.z, (int)marker.y, (int)marker.x);
    if ((double)intensity > (double)std::numeric_limits<int>::max())
        MCP3D_MESSAGE("warning: intensity value " + std::to_string(intensity) + " is greater than representable value of int type")
    marker_intensity = (int)intensity;
}

template <typename VType>
void App2ThresholdVolume(mcp3d::MImageBase &image, mcp3d::App2Parameters &app2_parameters)
{
    MCP3D_ASSERT(image.n_volumes() == 1)
    MCP3D_ASSERT(app2_parameters.foreground_percent() < 0 && app2_parameters.threshold_value() >= 0)
    int channel_number = image.loaded_view().view_channels()[0];
    if (sizeof(int) > sizeof(VType) &&
        (double)app2_parameters.threshold_value() > (double)std::numeric_limits<VType>::max())
        MCP3D_RUNTIME_ERROR("threshold value greater than representable value of the voxel type")
    mcp3d::ThresholdVolume<VType>(image.Volume<VType>(channel_number),
                                  image.loaded_view().xyz_dims(channel_number, 0),
                                  (VType)app2_parameters.threshold_value(), true);
}

template <typename VType>
void App2ThresholdVolumeTopPercentile(mcp3d::MImageBase &image,
                                      mcp3d::App2Parameters &app2_parameters)
{
    MCP3D_ASSERT(image.n_volumes() == 1)
    MCP3D_ASSERT(app2_parameters.foreground_percent() >= 0 && app2_parameters.foreground_percent() <= 1)
    int channel_number = image.loaded_view().view_channels()[0];
    VType threshold_value =
            mcp3d::ThresholdVolumeTopPercentile<VType>(image.Volume<VType>(channel_number),
                                                       image.loaded_view().xyz_dims(channel_number, 0),
                                                       app2_parameters.foreground_percent(), true);
    app2_parameters.set_threshold_value((int)threshold_value);
}

template <typename VType>
void App2ThresholdPlanesTopPercentile(mcp3d::MImageBase &image,
                                      mcp3d::App2Parameters &app2_parameters)
{
    MCP3D_ASSERT(image.n_volumes() == 1)
    MCP3D_ASSERT(app2_parameters.foreground_percent() >= 0 && app2_parameters.foreground_percent() <= 1)
    MCP3D_ASSERT(app2_parameters.foreground_percent_by_plane())
    int channel_number = image.loaded_view().view_channels()[0];
    std::vector<VType> plane_threshold_values =
            mcp3d::ThresholdPlanesTopPercentile<VType>(image.Volume<VType>(channel_number),
                                                       image.loaded_view().xyz_dims(channel_number, 0),
                                                       app2_parameters.foreground_percent(), 0, true);
    app2_parameters.set_plane_threshhold_values(std::vector<int>(plane_threshold_values.begin(),
                                                                 plane_threshold_values.end()));
}

#endif //MCP3D_APP2_IMAGE_PREPROCESS_HPP
