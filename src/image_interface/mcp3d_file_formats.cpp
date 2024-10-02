//
// Created by muyezhu on 3/17/19.
//
#include <string>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include "common/mcp3d_utility.hpp"
#include "mcp3d_file_formats.hpp"

using namespace std;

mcp3d::FileFormat mcp3d::StringToFileFormat(const string &format_str)
{
    if (format_str == "tif")
        return mcp3d::FileFormat::TIFF;
    else if (format_str == "ome.tif")
        return mcp3d::FileFormat::OMETIFF;
    else if (format_str == "ims")
        return mcp3d::FileFormat::IMARIS;
    else if (format_str == "hdf5")
        return mcp3d::FileFormat::HDF5;
    else
        return mcp3d::FileFormat::UNKNOWN;
}

string mcp3d::FileFormatToString(FileFormat format)
{
    if (format == mcp3d::FileFormat::TIFF)
        return "tif";
    else if (format == mcp3d::FileFormat::OMETIFF)
        return "ome.tif";
    else if (format == mcp3d::FileFormat::HDF5)
        return "hdf5";
    else if (format == mcp3d::FileFormat::IMARIS)
        return "ims";
    else
        return "unknown";
}

bool mcp3d::DirContainsFileFormat(const string& dir, mcp3d::FileFormat format)
{
    if (!mcp3d::IsDir(dir))
        MCP3D_OS_ERROR(dir + " is not a directory")
    if (format == mcp3d::FileFormat::IMARIS)
    {
        size_t n = mcp3d::StitchedImarisPathsInDir(dir).size();
        if (n == 1)
            return true;
        else if (n >= 2)
            MCP3D_RUNTIME_ERROR("more than one stitched imaris file found in " + dir)
    }
    else
    {
        boost::filesystem::directory_iterator dir_iter{dir};
        for(auto it = boost::filesystem::begin(dir_iter); it != boost::filesystem::end(dir_iter); ++it)
        {
            if (boost::filesystem::is_regular_file(it->path()))
            {
                string file_name = it->path().filename().string();
                if (mcp3d::FileIsFormat(file_name, format))
                    return true;
            }
        }
    }
    return false;
}

bool mcp3d::FileIsFormat(const string &file_path, FileFormat format)
{
    if (format == mcp3d::FileFormat::TIFF)
        return mcp3d::FileIsTiff(file_path);
    if (format == mcp3d::FileFormat::OMETIFF)
        return mcp3d::FileIsOmeTiff(file_path);
    if (format == mcp3d::FileFormat::IMARIS)
        return mcp3d::FileIsStitchedImaris(file_path);
    if (format == mcp3d::FileFormat::HDF5)
        return mcp3d::FileIsHdf5(file_path);
    if (format == mcp3d::FileFormat::UNKNOWN)
        return !mcp3d::FileIsTiff(file_path) && !mcp3d::FileIsOmeTiff(file_path) && !mcp3d::FileIsStitchedImaris(file_path) && !mcp3d::FileIsHdf5(file_path);
}

bool mcp3d::FileIsTiff(const string &file_path)
{
    if (mcp3d::FileIsOmeTiff(file_path))
        return false;
    return boost::algorithm::ends_with(file_path, mcp3d::TiffExtension());
}

bool mcp3d::FileIsOmeTiff(const std::string &file_path)
{
    return boost::algorithm::ends_with(file_path, mcp3d::OmeTiffExtension());
}

bool mcp3d::FileIsImaris(const string &file_path)
{
    return boost::algorithm::ends_with(file_path, mcp3d::ImarisExtension());
}

bool mcp3d::FileIsStitchedImaris(const string &file_path)
{
    if (!mcp3d::FileIsImaris(file_path))
        return false;
    return file_path.find("FusionStitcher") != string::npos;
}

bool mcp3d::FileIsHdf5(const string &file_path)
{
    return boost::algorithm::ends_with(file_path, mcp3d::Hdf5Extension());
}

vector<string> mcp3d::StitchedImarisPathsInDir(const std::string &dir)
{
    vector<string> stitched_ims_paths;
    if (!mcp3d::IsDir(dir))
    {
        MCP3D_MESSAGE(dir + " is not a directory")
        return stitched_ims_paths;
    }

    // find valid imaris file
    vector<string> file_paths = mcp3d::FilesInDir(dir);
    for (const string& file_path: file_paths)
        if (mcp3d::FileIsStitchedImaris(file_path))
            stitched_ims_paths.push_back(file_path);
    return stitched_ims_paths;
}

string mcp3d::StitchedImarisPathInDir(const std::string &dir)
{
    if (!mcp3d::IsDir(dir))
        MCP3D_RUNTIME_ERROR(dir + " is not a directory")
    vector<string> stitched_ims_paths(mcp3d::StitchedImarisPathsInDir(dir));
    if (stitched_ims_paths.empty())
        MCP3D_RUNTIME_ERROR("can not find stitched imaris file under " + dir + ": file name must contain FusionStitcher and ends with ims")
    if (stitched_ims_paths.size() > 1)
        MCP3D_RUNTIME_ERROR("multiple stitched imaris files under " + dir)
    return stitched_ims_paths[0];
}