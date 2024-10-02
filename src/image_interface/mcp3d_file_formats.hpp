//
// Created by muyezhu on 3/17/19.
//

#ifndef MCP3D_MCP3D_FORMATS_HPP
#define MCP3D_MCP3D_FORMATS_HPP

#include <vector>
#include <unordered_set>


namespace mcp3d
{

enum class FileFormat
{
    TIFF, OMETIFF, IMARIS, HDF5, UNKNOWN = -1
};

inline std::string ImarisExtension()  { return ".ims"; }

inline std::string Hdf5Extension()  { return ".hdf5"; }

inline std::string TiffExtension()  { return ".tif"; }

inline std::string OmeTiffExtension()  { return "ome.tif"; }

FileFormat StringToFileFormat(const std::string &format_str);

std::string FileFormatToString(FileFormat format);

/// does not assert file_path exists
/// IMARIS format is true when FileIsStitchedImaris returns true (instead of FileIsImaris)
/// if format is UNKNOWN, return true if file_path is of none of the known format
bool FileIsFormat(const std::string& file_path, FileFormat format);

/// does not assert file_path exists
bool FileIsTiff(const std::string &file_path);

/// does not assert file_path exists
bool FileIsOmeTiff(const std::string &file_path);

/// does not assert file_path exists
bool FileIsImaris(const std::string& file_path);

/// does not assert file_path exists. return true if file path ends with .ims and contains FusionStitcher
bool FileIsStitchedImaris(const std::string& file_path);

/// does not assert file_path exists
bool FileIsHdf5(const std::string &file_path);

std::vector<std::string> StitchedImarisPathsInDir(const std::string &dir);

/// throws error if dir is not a directory, or if dir contains none or more than
/// one stitched imaris files
std::string StitchedImarisPathInDir(const std::string &dir);

/// use a directory iterator and search for files consistent with format. return true whenever
/// a consistent file is found. for imaris files, throw exception if more than one stitched imaris files are found.
/// return false if no consistent file is found
bool DirContainsFileFormat(const std::string& dir, FileFormat format);

inline std::vector<FileFormat> FileFormatSearchOrder()
{ return { FileFormat::IMARIS, FileFormat::HDF5, FileFormat::TIFF, FileFormat::OMETIFF }; }

}


#endif //MCP3D_MCP3D_FORMATS_HPP
