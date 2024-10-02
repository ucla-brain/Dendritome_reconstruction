//
// Created by muyezhu on 5/24/19.
//

#ifndef MCP3D_VAA3D_APP2_HPP
#define MCP3D_VAA3D_APP2_HPP

#include <iostream>
#include <sstream>
#include "common/mcp3d_utility.hpp"

namespace mcp3d
{

class App2Parameters
{
public:
    App2Parameters(): foreground_percent_by_plane_(false), gsdt_(false), 
                      coverage_prune_(true), allow_gap_(false), radius_from_2d_(true),
                      swc_resample_(true), high_intensity_(false), brightfield_(false),
                      threshold_value_(-1), cnn_type_(2), plane_threshhold_values_(std::vector<int>{}),
                      foreground_percent_(-0.01), sr_ratio_(1.0 / 3),
                      length_thresh_(5.0), marker_file_path_(std::string{}) {}

    std::string MetaString();
    
    // getters
    bool foreground_percent_by_plane() const  { return foreground_percent_by_plane_; }
    bool gsdt() const { return gsdt_; }
    bool coverage_prune() const  { return coverage_prune_; }
    bool allow_gap() const  { return allow_gap_; }
    int threshold_value() const  { return threshold_value_; }
    std::vector<int> plane_threshold_values() const  { return plane_threshhold_values_; }
    std::string plane_threshold_values_string() const
    { return mcp3d::JoinVector(plane_threshhold_values_, ", ", true); }
    double foreground_percent() const  { return foreground_percent_; }
    double length_thresh() const  { return length_thresh_; }
    int cnn_type() const  { return cnn_type_; }
    double sr_ratio() const  { return sr_ratio_; }
    bool cube_256() const  { return cube_256_; }
    bool radius_from_2d() const  { return radius_from_2d_; }
    bool swc_resample() const  { return swc_resample_; }
    bool high_intensity() const  { return high_intensity_; }
    bool brightfield() const  { return brightfield_; }
    std::string marker_file_path() const  { return marker_file_path_; }
    // setters
    void set_foreground_percent_by_plane(bool foreground_percent_by_plane)  { foreground_percent_by_plane_ = foreground_percent_by_plane; }
    void set_gsdt(bool is_gsdt) { gsdt_ = is_gsdt; }
    void seg_coverage_prune(bool coverage_prune)  { coverage_prune_ = coverage_prune; }
    void set_allow_gap(bool allow_gap)  { allow_gap_ = allow_gap; }
    void set_threshold_value(int threshold_value)  { threshold_value_ = std::max(threshold_value, 0); }
    void set_plane_threshhold_values(const std::vector<int>& plane_threshhold_values)  { plane_threshhold_values_ = plane_threshhold_values; }
    void set_foreground_percent(double foreground_percent)  { foreground_percent_ = std::max(0.0, std::min(1.0, foreground_percent)); }
    void set_length_thresh(double length_thresh)  { length_thresh_ = std::max(length_thresh, 1.0); }
    void set_cnn_type(int cnn_type)  { if (cnn_type == 1 || cnn_type == 2 || cnn_type == 3) cnn_type_ = cnn_type; }
    void set_sr_ratio(double sr_ratio)  { sr_ratio_ = std::max(1e-9, std::min(1.0, sr_ratio)); }
    void set_cube_256(bool cube_256)  { cube_256_ = cube_256; }
    void set_radius_from_2d(bool radius_from_2d)  { radius_from_2d_ = radius_from_2d; }
    void set_swc_resample(bool swc_resample)  { swc_resample_ = swc_resample; }
    void set_high_intensity(bool high_intensity)  { high_intensity_ = high_intensity; }
    void set_brightfield(bool brightfield)  { brightfield_ = brightfield; }
    void set_marker_file_path(const std::string& marker_file_path)
    {
        if (mcp3d::IsFile(marker_file_path))
            marker_file_path_ = marker_file_path;
        else
        MCP3D_MESSAGE(marker_file_path + " is not a file")
    }

private:
    bool foreground_percent_by_plane_, gsdt_, coverage_prune_, allow_gap_,
         cube_256_, radius_from_2d_, swc_resample_,
         high_intensity_, brightfield_;
    int threshold_value_, cnn_type_;
    std::vector<int> plane_threshhold_values_;
    double foreground_percent_, sr_ratio_, length_thresh_;
    std::string marker_file_path_;
};

class App2CommandLineArgs
{
public:
    App2CommandLineArgs(): app2_parameters_(App2Parameters{}),
                           image_root_dir_(std::string()), swc_path_(std::string()),
                           channel_number_(0), resolution_level_(0),
                           image_offsets_({0, 0, 0}), image_extents_({0, 0, 0}) {}
    static void PrintUsage();
    std::string MetaString();
    void PrintParameters()  { std::cout << MetaString() << std::endl; }

    App2Parameters& app2_parameters()  { return app2_parameters_; }
    // getters
    const App2Parameters& app2_parameters() const  { return app2_parameters_; }
    std::string image_root_dir() const  { return image_root_dir_; }
    std::string swc_path() const  { return swc_path_; }
    int channel_number() const  { return channel_number_; }
    int resolution_level() const  { return resolution_level_; }
    std::vector<int> image_offsets() const { return image_offsets_; }
    std::vector<int> image_extents() const { return image_extents_; }

    // setters
    void set_image_root_dir(const std::string& image_root_dir)
    {
        if (!mcp3d::IsDir(image_root_dir))
            MCP3D_INVALID_ARGUMENT(image_root_dir + " is not a directory")
        image_root_dir_ = image_root_dir;
    }

    void set_swc_path(const std::string& swc_path) { swc_path_ = swc_path; }

    void set_channel_number(int channel_number)
    {
        if (channel_number < 0)
            MCP3D_INVALID_ARGUMENT("channel number must be non negative, received " + std::to_string(channel_number))
        channel_number_ = channel_number;
    }

    void set_resolution_level(int resolution_level)
    {
        if (resolution_level < 0)
            MCP3D_INVALID_ARGUMENT("resolution_level must be non negative, received " + std::to_string(resolution_level))
        resolution_level_ = resolution_level;
    }

    void set_image_offsets(const std::vector<int>& image_offsets)
    {
        for (size_t i = 0; i < std::min((size_t)3, image_offsets.size()); ++i)
        {
            if (image_offsets[i] < 0)
                MCP3D_INVALID_ARGUMENT("offsets values must be non negative")
            image_offsets_[i] = image_offsets[i];
        }
    }

    void set_image_extents(const std::vector<int>& image_extents)
    {
        for (size_t i = 0; i < std::min((size_t)3, image_extents.size()); ++i)
        {
            if (image_extents[i] < 0)
                MCP3D_INVALID_ARGUMENT("offsets values must be non negative")
            image_extents_[i] = image_extents[i];
        }
    }

private:
    App2Parameters app2_parameters_;
    std::string image_root_dir_, swc_path_;
    int channel_number_, resolution_level_;
    std::vector<int> image_offsets_, image_extents_;
};

bool ParseApp2Args(int argc, char * argv[], App2CommandLineArgs& args);

}


#endif //MCP3D_VAA3D_APP2_HPP
