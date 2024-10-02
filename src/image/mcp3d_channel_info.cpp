//
// Created by muyezhu on 2/11/18.
//
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include "image_layout/mcp3d_image_layout_constants.hpp"
#include "mcp3d_image_constants.hpp"
#include "mcp3d_image_utils.hpp"
#include "mcp3d_channel_info.hpp"

using namespace std;
using json = nlohmann::json;


mcp3d::MChannelPyrInfo::MChannelPyrInfo(const json &info_json, const shared_ptr<mcp3d::MChannelPyrSlices>& channel_pyr_slices):
                                     channel_pyr_slices_(channel_pyr_slices)
{
    try
    {
        ParseJson(info_json);
        AssertValid();
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

mcp3d::MChannelPyrInfo::MChannelPyrInfo(const mcp3d::MChannelPyrInfo& other)
{
    channel_pyr_slices_ = make_shared<mcp3d::MChannelPyrSlices>(other.channel_pyr_slices_->channel_pyr_dir());
    file_format_ = other.file_format_;
    resolution_level_ = other.resolution_level_;
    xdim_ = other.xdim_;
    ydim_ = other.ydim_;
    zdim_ = other.zdim_;
    chunk_xdim_ = other.chunk_xdim_;
    chunk_ydim_ = other.chunk_ydim_;
    chunk_zdim_ = other.chunk_zdim_;
    channel_pyr_dir_ = other.channel_pyr_dir_;
    slice_image_names_ = other.slice_image_names_;
    voxel_type_ = other.voxel_type_;
}

mcp3d::MChannelPyrInfo::MChannelPyrInfo(mcp3d::MChannelPyrInfo&& other) noexcept
{
    // move assignment
    channel_pyr_slices_ = other.channel_pyr_slices_;
    file_format_ = other.file_format_;
    resolution_level_ = other.resolution_level_;
    xdim_ = other.xdim_;
    ydim_ = other.ydim_;
    zdim_ = other.zdim_;
    chunk_xdim_ = other.chunk_xdim_;
    chunk_ydim_ = other.chunk_ydim_;
    chunk_zdim_ = other.chunk_zdim_;
    channel_pyr_dir_ = other.channel_pyr_dir_;
    slice_image_names_ = other.slice_image_names_;
    voxel_type_ = other.voxel_type_;
}

bool mcp3d::MChannelPyrInfo::operator==(const mcp3d::MChannelPyrInfo &other) const
{
    return  *(channel_pyr_slices_.get()) == *(other.channel_pyr_slices_.get()) &&
            file_format_ == other.file_format_ &&
            resolution_level_ == other.resolution_level_ &&
            xdim_ == other.xdim_ &&
            ydim_ == other.ydim_ &&
            zdim_ == other.zdim_ &&
            chunk_xdim_ == other.chunk_xdim_ &&
            chunk_ydim_ == other.chunk_ydim_ &&
            chunk_zdim_ == other.chunk_zdim_ &&
            channel_pyr_dir_ == other.channel_pyr_dir_ &&
            slice_image_names_ == other.slice_image_names_ &&
            voxel_type_ == other.voxel_type_;
}

void mcp3d::MChannelPyrInfo::Clear()
{
    channel_pyr_slices_.reset();
    file_format_ = mcp3d::FileFormat::UNKNOWN;
    resolution_level_ = -1;
    xdim_ = 0;
    ydim_ = 0;
    zdim_ = 0;
    chunk_xdim_ = 0;
    chunk_ydim_ = 0;
    chunk_zdim_ = 0;
    channel_pyr_dir_.clear();
    slice_image_names_.clear();
    voxel_type_ = mcp3d::VoxelType::UNKNOWN;
    volume_complete_ = mcp3d::VolumeComplete::UNKNOWN;
}

bool mcp3d::MChannelPyrInfo::empty() const
{
    if (!channel_pyr_slices_)
        return true;
    return channel_pyr_slices_->empty();
}

void mcp3d::MChannelPyrInfo::UpdateVolumeCompleteness()
{
    if (!DimensionsDefined())
        return;
    for (int i = 0; i < slice_number(mcp3d::ChannelAxis::Z); ++i)
        for (int j = 0; j < slice_number(mcp3d::ChannelAxis::Y); ++j)
            for (int k = 0; k < slice_number(mcp3d::ChannelAxis::X); ++k)
                if (NumberOfImageChunksInSlice(i, j, k) < ExpectedNumberOfImageChunksInSlice(i, j, k))
                {
                    volume_complete_ = VolumeComplete::INCOMPLETE;
                    return;
                }
    volume_complete_ = VolumeComplete::COMPLETE;
}

void mcp3d::MChannelPyrInfo::ReadSliceImageNames()
{
    if (empty())
        return;
    slice_image_names_.clear();
    for (const auto& slice_name: channel_pyr_slices_->SliceNames())
    {
        slice_image_names_.insert(make_pair(slice_name, vector<string>{}));
        if (file_format_ == mcp3d::FileFormat::UNKNOWN)
            continue;
        string slice_dir(mcp3d::JoinPath(channel_pyr_dir_, slice_name));
        boost::filesystem::path p(slice_dir);
        boost::filesystem::directory_iterator it_end;
        for(auto it = boost::filesystem::directory_iterator(p); it != it_end; ++it)
        {
            if (boost::filesystem::is_regular_file(it->status()))
            {
                string file_name = it->path().filename().string();
                if (mcp3d::FileIsFormat(file_name, file_format_))
                    slice_image_names_.at(slice_name).push_back(file_name);
            }
        }
        sort(slice_image_names_.at(slice_name).begin(), slice_image_names_.at(slice_name).end());
    }
    UpdateVolumeCompleteness();
}

mcp3d::FileFormat mcp3d::MChannelPyrInfo::ParseFileFormat(const json &info_json)
{
    try
    {
        mcp3d::FileFormat parsed_file_format = mcp3d::StringToFileFormat(info_json["format"]);
        return parsed_file_format;
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR("file format parsing failure\n" + string(e.what()))
    }
}

int mcp3d::MChannelPyrInfo::ParseResolutionLevel(const nlohmann::json &info_json)
{
    try
    {
        int parsed_resolution_level = info_json["resolution level"];
        return parsed_resolution_level;
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR("resolution level parsing failure\n" + string(e.what()))
    }
}

string mcp3d::MChannelPyrInfo::ParseChannelPyrDir(const json& info_json)
{
    try
    {
        string parsed_channel_pyr_dir = info_json["channel pyramid directory"];
        return mcp3d::RemoveTrailingSeparator(parsed_channel_pyr_dir);
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR("channel pyramid directory parsing failure\n" + string(e.what()))
    }
}

int mcp3d::MChannelPyrInfo::ParseDim(const json& info_json, mcp3d::ChannelAxis axis)
{
    string axis_str = mcp3d::ChannelAxisStr(axis);
    axis_str.append(" dimension");
    try
    {
        int parsed_dim = info_json[axis_str];
        return parsed_dim;
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR(axis_str + " parsing failure\n" + string(e.what()))
    }
}

int mcp3d::MChannelPyrInfo::ParseChunkDim(const json &info_json, mcp3d::ChannelAxis axis)
{
    string key = "chunk ";
    key.append(mcp3d::ChannelAxisStr(axis)).append(" dimension");
    try
    {
        int parsed_chunk_dim = info_json[key];
        return parsed_chunk_dim;
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR(key + " parsing failure\n" + string(e.what()))
    }
}

mcp3d::VoxelType mcp3d::MChannelPyrInfo::ParseVoxelType(const json &info_json)
{
    try
    {
        string parsed_voxel_type = info_json["voxel type"];
        return mcp3d::StringToVoxelType(parsed_voxel_type);
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR("voxel type parsing failure\n" + string(e.what()))
    }
}

void mcp3d::MChannelPyrInfo::ParseSliceImageNames(const json &info_json)
{
    try
    {
        slice_image_names_ = info_json["image sequence"].get<unordered_map<string, vector<string>>>();
    }
    catch (const exception& e)
    {
        MCP3D_RUNTIME_ERROR("image sequence parsing failure\n" + string(e.what()))
    }
}

void mcp3d::MChannelPyrInfo::ParseJson(const json& info_json)
{
    // image format
    file_format_ = ParseFileFormat(info_json);
    // resolution level
    resolution_level_ = ParseResolutionLevel(info_json);
    // channel pyramid level directory
    channel_pyr_dir_ = ParseChannelPyrDir(info_json);

    // bool operator of shared_ptr
    if (!channel_pyr_slices_)
        channel_pyr_slices_ = make_shared<mcp3d::MChannelPyrSlices>(channel_pyr_dir_);

    // image dimensions
    xdim_ = ParseDim(info_json, mcp3d::ChannelAxis::X);
    ydim_ = ParseDim(info_json, mcp3d::ChannelAxis::Y);
    zdim_ = ParseDim(info_json, mcp3d::ChannelAxis::Z);

    // chunk dimensions
    chunk_xdim_ = ParseChunkDim(info_json, mcp3d::ChannelAxis::X);
    chunk_ydim_ = ParseChunkDim(info_json, mcp3d::ChannelAxis::Y);
    chunk_zdim_ = ParseChunkDim(info_json, mcp3d::ChannelAxis::Z);

    // voxel data type
    voxel_type_ = ParseVoxelType(info_json);

    // image name sequence in folder slices
    ParseSliceImageNames(info_json);
}

void mcp3d::MChannelPyrInfo::AssertValid() const
{
    if (empty())
        return;
    MCP3D_ASSERT(voxel_type_ != mcp3d::VoxelType::UNKNOWN)
    MCP3D_ASSERT(channel_pyr_dir_ == channel_pyr_slices_->channel_pyr_dir())
    MCP3D_ASSERT(mcp3d::IsDir(channel_pyr_dir_))
    MCP3D_ASSERT(file_format_ == channel_pyr_slices_->file_format())
    if (file_format_ != mcp3d::FileFormat::IMARIS)
        MCP3D_ASSERT(resolution_level_ == mcp3d::MChannelLayout::DirPyrLevel(channel_pyr_dir_))
    vector<string> slice_names{channel_pyr_slices_->SliceNames().size()};
    MCP3D_ASSERT(slice_image_names_.size() == slice_names.size())
    for (const auto& slice_name: slice_names)
        MCP3D_ASSERT(slice_image_names_.find(slice_name) != slice_image_names_.end())
    if (!DimensionsDefined())
        return;
    for (const auto& axis: vector<mcp3d::ChannelAxis>({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X}))
    {
        MCP3D_ASSERT(dim(axis) > 0 && chunk_dim(axis) > 0 && slice_dim(axis) > 0)
        MCP3D_ASSERT(dim(axis) >= slice_dim(axis) && dim(axis) >= chunk_dim(axis))
        MCP3D_ASSERT(slice_dim(axis) >= chunk_dim(axis) && slice_dim(axis) % chunk_dim(axis) == 0)
        MCP3D_ASSERT(chunk_dim(axis) * n_chunks(axis) >= dim(axis) && chunk_dim(axis) * (n_chunks(axis) - 1) < dim(axis))
        MCP3D_ASSERT(slice_dim(axis) * slice_dim(axis) >= chunk_dim(axis) * n_chunks(axis))
    }
    MCP3D_ASSERT(volume_complete_ != mcp3d::VolumeComplete::UNKNOWN)
    for (int i = 0; i < slice_number(mcp3d::ChannelAxis::Z); ++i)
        for (int j = 0; j < slice_number(mcp3d::ChannelAxis::Z); ++j)
            for (int k = 0; k < slice_number(mcp3d::ChannelAxis::Z); ++k)
            {
                if (volume_complete_ == mcp3d::VolumeComplete::COMPLETE)
                    MCP3D_ASSERT(NumberOfImageChunksInSlice(i, j, k) == ExpectedNumberOfImageChunksInSlice(i, j, k))
                else
                    MCP3D_ASSERT(NumberOfImageChunksInSlice(i, j, k) <= ExpectedNumberOfImageChunksInSlice(i, j, k))
            }
}

json mcp3d::MChannelPyrInfo::Jsonize() const
{
    json info_json;
    info_json["format"] = mcp3d::FileFormatToString(file_format_);
    info_json["resolution level"] = resolution_level_;
    info_json["channel pyramid directory"] = channel_pyr_dir_;
    info_json["x dimension"] = xdim_;
    info_json["y dimension"] = ydim_;
    info_json["z dimension"] = zdim_;
    info_json["chunk x dimension"] = chunk_xdim_;
    info_json["chunk y dimension"] = chunk_ydim_;
    info_json["chunk z dimension"] = chunk_zdim_;
    info_json["image sequence"] = slice_image_names_;
    info_json["voxel type"] = mcp3d::VoxelTypeToString(voxel_type_);
    return info_json;
}

int mcp3d::MChannelPyrInfo::NumberOfImageChunksInSlice(int zslice_id, int yslice_id, int xslice_id) const
{
    string slice_name = channel_pyr_slices_->SliceNameFromSliceIDs(zslice_id, yslice_id, xslice_id);
    return (int)slice_image_names_.at(slice_name).size();
}

int mcp3d::MChannelPyrInfo::ExpectedNumberOfImageChunksInSlice(int zslice_id, int yslice_id, int xslice_id) const
{
    return ExpectedNumberOfImageChunksInSliceAlongAxis(ChannelAxis::Z, zslice_id) *
           ExpectedNumberOfImageChunksInSliceAlongAxis(ChannelAxis::Y, yslice_id) *
           ExpectedNumberOfImageChunksInSliceAlongAxis(ChannelAxis::X, xslice_id);
}

int mcp3d::MChannelPyrInfo::ExpectedNumberOfImageChunksInSliceAlongAxis(ChannelAxis axis, int slice_id) const
{
    if (!DimensionsDefined())
        return 0;
    MCP3D_ASSERT(slice_id >= 0 && slice_id < slice_number(axis))
    // if the axis is covered by a single slice, return 1
    if (slice_number(axis) == 1)
        return 1;
    if (slice_id < slice_number(axis) - 1)
        return slice_dim(axis) / chunk_dim(axis);
    else
    {
        int dim_residual = dim(axis) % slice_dim(axis);
        return dim_residual / chunk_dim(axis) + (int)(dim_residual % chunk_dim(axis) == 0);
    }
}

int mcp3d::MChannelPyrInfo::TotalNumberOfImageChunks() const
{
    int n = 0;
    for (const auto& item: slice_image_names_)
        n += (int)item.second.size();
    return n;
}

int mcp3d::MChannelPyrInfo::slice_dim(ChannelAxis axis) const
{
    if (channel_pyr_slices_->slice_dim(axis) == mcp3d::MAX_AXIS_SLICE_DIM)
        return dim(axis);
    return channel_pyr_slices_->slice_dim(axis);
}

string mcp3d::MChannelPyrInfo::ImagePath(int z, int y, int x) const
{
    MCP3D_ASSERT(DimensionsDefined())
    MCP3D_ASSERT(volume_complete_ == mcp3d::VolumeComplete::COMPLETE)
    if (!ZyxInVolume(z, y, x))
        return string{};
    string slice_name(channel_pyr_slices_->SliceNameFromCoordinates(z, y, x));
    int slice_zid = z / slice_dim(mcp3d::ChannelAxis::Z),
        slice_yid = y / slice_dim(mcp3d::ChannelAxis::Y),
        slice_xid = x / slice_dim(mcp3d::ChannelAxis::X);
    int slice_zval = z % slice_dim(mcp3d::ChannelAxis::Z),
        slice_yval = y % slice_dim(mcp3d::ChannelAxis::Y),
        slice_xval = x % slice_dim(mcp3d::ChannelAxis::X);
    int chunk_zid = slice_zval / chunk_zdim_, chunk_yid = slice_yval / chunk_ydim_, chunk_xid = slice_xval / chunk_xdim_;
    long chunk_index = mcp3d::LinearAddress(vector<int>({ExpectedNumberOfImageChunksInSliceAlongAxis(mcp3d::ChannelAxis::Z, slice_zid),
                                                         ExpectedNumberOfImageChunksInSliceAlongAxis(mcp3d::ChannelAxis::Y, slice_yid),
                                                         ExpectedNumberOfImageChunksInSliceAlongAxis(mcp3d::ChannelAxis::X, slice_xid)}),
                                            chunk_zid, chunk_yid, chunk_xid);
    return mcp3d::JoinPath(channel_pyr_dir_, slice_name, slice_image_names_.at(slice_name)[chunk_index]);
}

int mcp3d::MChannelPyrInfo::dim(ChannelAxis axis) const
{
    if (axis == mcp3d::ChannelAxis::Z)
        return zdim();
    else if (axis == mcp3d::ChannelAxis::Y)
        return ydim();
    else
        return xdim();
}

int mcp3d::MChannelPyrInfo::chunk_dim(ChannelAxis axis) const
{
    if (axis == mcp3d::ChannelAxis::Z)
        return chunk_zdim();
    else if (axis == mcp3d::ChannelAxis::Y)
        return chunk_ydim();
    else
        return chunk_xdim();
}

int mcp3d::MChannelPyrInfo::n_chunks(ChannelAxis axis) const
{
    if (axis == mcp3d::ChannelAxis::Z)
        return n_zchunks();
    else if (axis == mcp3d::ChannelAxis::Y)
        return n_ychunks();
    else
        return n_xchunks();
}

mcp3d::MChannelInfo::MChannelInfo(const mcp3d::MChannelInfo& other): channel_pyr_infos_(std::unordered_map<int, MChannelPyrInfo> {})
{
    channel_layout_ = make_shared<mcp3d::MChannelLayout>(other.channel_layout());
    for (const auto& item: other.channel_pyr_infos_)
    {
        int pyr_level = item.first;
        AddPyrInfo(mcp3d::MChannelPyrInfo(other.channel_pyr_infos_.at(pyr_level).Jsonize(), channel_layout_->channel_pyr_slices_ptr(pyr_level)));
        channel_pyr_infos_[pyr_level].AssertValid();
    }
}

void mcp3d::MChannelInfo::AddPyrInfo(MChannelPyrInfo &&pyr_info)
{
    if (pyr_info.empty())
        return;
    channel_layout_->ReadChannelLayout();
    int pyr_level = pyr_info.resolution_level();
    MCP3D_ASSERT(channel_layout_->HasPyrLevel(pyr_level))
    // assert equality of MChannelPyrSlices managed by channel_layout_ and by pyr_info
    MCP3D_ASSERT(channel_layout_->channel_pyr_slices(pyr_level) == pyr_info.channel_pyr_slices())
    channel_pyr_infos_.emplace(pyr_level, move(pyr_info));
    // copy assign shared_ptr<MChannelPyrSlices> under channel_layout_ to channel_pyr_infos_[pyr_level]
    channel_pyr_infos_.at(pyr_level).channel_pyr_slices_ = channel_layout_->channel_pyr_slices_ptr(pyr_level);
}

void mcp3d::MChannelInfo::RefreshChannelLayout()
{
    channel_layout_->ReadChannelLayout();
    for (auto& item: channel_pyr_infos_)
    {
        if (!channel_layout_->HasPyrLevel(item.first))
            RemovePyrInfo(item.first);
    }
}

bool mcp3d::MChannelInfo::LevelPathsExist(int pyr_level) const
{
    MCP3D_ASSERT(pyr_level >= 0 && pyr_level < n_pyr_infos())
    bool is_mpi = mcp3d::MPIInitialized();
    int n_threads = is_mpi? 1 : mcp3d::DefaultNumThreads();
    const MChannelPyrInfo& pyr_info = channel_pyr_infos().at(pyr_level);
    int n_slices_per_thread = TotalSliceNumber(pyr_level) / n_threads;
    bool exist = true;
    vector<string> slice_names = channel_pyr_slices(pyr_level).SliceNames();
#pragma omp parallel num_threads(n_threads)
    {
        try
        {
            CHECK_PARALLEL_MODEL
            int thread_id = omp_get_thread_num();
            int slice_id_begin = thread_id * n_slices_per_thread,
                    slice_id_end = min(slice_id_begin + n_slices_per_thread,
                                       TotalSliceNumber(pyr_level));
            for (int j = slice_id_begin; j < slice_id_end; ++j)
            {
                string slice_name = slice_names[j];
                string slice_dir = mcp3d::JoinPath(pyr_info.channel_pyr_dir_, slice_name);
                for (const auto& image_name: pyr_info.slice_image_names_.at(slice_name))
                    if (!mcp3d::IsFile(mcp3d::JoinPath(slice_dir, image_name)))
                        MCP3D_RUNTIME_ERROR(mcp3d::JoinPath(slice_dir, image_name) + " is not a file")
            }
        }
        catch (const mcp3d::MCPRuntimeError& e)
        {
            mcp3d::PrintNested(e);
            exist = false;
        }
    }
    return exist;
}

void mcp3d::MChannelInfo::Save() const
{
    if (n_pyr_infos() == 0)
    {
        MCP3D_MESSAGE("no MChannelPyrInfo associated with MChannelInfo instance. do nothing.")
        return;
    }
    json pyr_infos;
    for (const auto pyr_level: pyr_info_levels())
        pyr_infos[mcp3d::MChannelLayout::PyrLevelDirName(pyr_level)] = channel_pyr_infos_.at(pyr_level).Jsonize();
    ofstream out_file(channel_info_path(), fstream::out);
    out_file << pyr_infos.dump(4) << endl;
    out_file.close();
}

void mcp3d::MChannelInfo::Load()
{
    // hold content of channel_pyr_infos_ in previous_channel_pyr_infos
    unordered_map<int, mcp3d::MChannelPyrInfo> previous_channel_pyr_infos;
    previous_channel_pyr_infos.swap(channel_pyr_infos_);
    // read and construct MChannelPyrInfo instances from __channel_info__.json
    if (!mcp3d::IsFile(channel_info_path()))
    {
        MCP3D_MESSAGE("__channel_info__.json not found under " + channel_root_dir() + ". do nothing.")
        return;
    }
    json pyr_info_jsons;
    ifstream ifs(channel_info_path());
    ifs >> pyr_info_jsons;
    ifs.close();
    json pyr_info_json;
    for (auto& item: pyr_info_jsons.items())
    {
        pyr_info_json = item.value();
        int pyr_level = mcp3d::MChannelPyrInfo::ParseResolutionLevel(pyr_info_json);
        AddPyrInfo(mcp3d::MChannelPyrInfo{pyr_info_json, channel_layout_->channel_pyr_slices_ptr(pyr_level)});
    }
    try
    {
        // if loaded MChannelPyrInfo instances consititute valid MChannelInfo instance,
        // clear previous_channel_pyr_info
        AssertValid();
    }
    catch (...)
    {
        // otherwise reinstate channel_pyr_infos_ to previous content
        MCP3D_PRINT_NESTED_EXCEPTION
        cout << "loading from __channel_info__.json will result in invalid instace. do nothing." << endl;
        channel_pyr_infos_.swap(previous_channel_pyr_infos);
    }
    previous_channel_pyr_infos.clear();
}

vector<int> mcp3d::MChannelInfo::ZyxTransform(int pyr_level, int z, int y, int x, bool global_to_local) const
{
    if (pyr_level == 0)
        return vector<int>({z, y, x});
    vector<int> local_zyx;
    unordered_map<int, int> ratios(pyr_ratios(mcp3d::ChannelAxis::Z));
    if (ratios.find(pyr_level) == ratios.end())
    {
        MCP3D_MESSAGE("can not translate global coordinates to local coordinates at pyr level " +
                      to_string(pyr_level) + ": required MChannelPyrInfo instance(s) not found.")
        return local_zyx;
    }
    local_zyx.push_back(global_to_local ? z / ratios.at(pyr_level) : z * ratios.at(pyr_level));
    ratios = pyr_ratios(mcp3d::ChannelAxis::Y);
    local_zyx.push_back(global_to_local ? y / ratios.at(pyr_level) : y * ratios.at(pyr_level));
    local_zyx.push_back(global_to_local ? x / ratios.at(pyr_level) : x * ratios.at(pyr_level));
    return local_zyx;
}

vector<int> mcp3d::MChannelInfo::pyr_info_levels() const
{
    vector<int> levels;
    for (const auto& item: channel_pyr_infos_)
        levels.push_back(item.first);
    sort(levels.begin(), levels.end());
    return levels;
}

const mcp3d::MChannelPyrInfo& mcp3d::MChannelInfo::channel_pyr_info(int pyr_level) const
{
    MCP3D_ASSERT(HasPyrInfo(pyr_level))
    return channel_pyr_infos_.at(pyr_level);
}

const mcp3d::MChannelPyrSlices& mcp3d::MChannelInfo::channel_pyr_slices(int pyr_level) const
{
    MCP3D_ASSERT(HasPyrLevel(pyr_level))
    return channel_layout_->channel_pyr_slices(pyr_level);
}

unordered_map<int, int> mcp3d::MChannelInfo::pyr_ratios(mcp3d::ChannelAxis axis) const
{
    unordered_map<int, int> ratios;
    if (!HasPyrInfo(0))
    {
        MCP3D_MESSAGE("no level 0 MChannelPyrInfo found. can not calculate pyramid ratios")
        return ratios;
    }
    const mcp3d::MChannelPyrInfo& pyr_info_0 = channel_pyr_infos_.at(0);
    for (const auto& level_pyr_info: channel_pyr_infos_)
    {
        int pyr_level = level_pyr_info.first;
        const mcp3d::MChannelPyrInfo& pyr_info = level_pyr_info.second;
        if (axis == mcp3d::ChannelAxis::X)
            ratios[pyr_level] = pyr_info_0.xdim() / pyr_info.xdim();
        else if (axis == mcp3d::ChannelAxis::Y)
            ratios[pyr_level] = pyr_info_0.ydim() / pyr_info.ydim();
        else
            ratios[pyr_level] = pyr_info_0.zdim() / pyr_info.zdim();
    }
    return ratios;
}

int mcp3d::MChannelInfo::z_scale_start_level() const
{
    int z = mcp3d::ZSCALE_NONE;
    unordered_map<int, int> ratios(pyr_ratios(mcp3d::ChannelAxis::Z));
    for (const auto& item: ratios)
    {
        if (item.second > 1)
            z = z < 0 ? item.first : min(z, item.first);
    }
    return z;
}

void mcp3d::MChannelInfo::AssertValid() const
{
    channel_layout_->AssertValid();
    // if no MImagePyrInfo instances exist, exit
    if (n_pyr_infos() == 0)
        return;
    // assert all levels in channel_pyr_infos_ exist in channel_layout_
    for (const auto pyr_level: pyr_info_levels())
    {
        MCP3D_ASSERT(HasPyrLevel(pyr_level))
        MCP3D_ASSERT(HasPyrInfo(pyr_level))
    }
    // assert MChannelPyrInfo instances are valid. channel_layout_ and channel_pyr_infos_ manages same MChannelPyrSlices objects
    for (const auto pyr_level: pyr_info_levels())
    {
        MCP3D_ASSERT(!channel_pyr_infos_.at(pyr_level).empty())
        channel_pyr_infos_.at(pyr_level).AssertValid();
        MCP3D_ASSERT(channel_layout_->channel_pyr_slices_ptr(pyr_level).get() == channel_pyr_infos_.at(pyr_level).channel_pyr_slices_ptr().get())
    }
    // validate pyramid ratios if level 0 MChannelPyrInfo exists
    if (!HasPyrLevel(0))
        return;
    // pyramid xy ratios correct
    unordered_map<int, int> xratios = pyr_ratios(mcp3d::ChannelAxis::X);
    for (const auto pyr_level: pyr_info_levels())
    {
        MCP3D_ASSERT (xratios.at(pyr_level) == mcp3d::IntPow<int>(2, pyr_level));
        MCP3D_ASSERT(xratios.at(pyr_level) <= xdim(pyr_level))
    }
    unordered_map<int, int> yratios = pyr_ratios(mcp3d::ChannelAxis::Y);
    for (const auto pyr_level: pyr_info_levels())
    {
        MCP3D_ASSERT (yratios.at(pyr_level) == mcp3d::IntPow<int>(2, pyr_level));
        MCP3D_ASSERT(yratios.at(pyr_level) <= ydim(pyr_level))
    }
    // pyramid z ratios correct. ratio[pyr_level] = ratio[pyr_level - 1] or ratio[pyr_level] * 2
    unordered_map<int, int> zratios = pyr_ratios(mcp3d::ChannelAxis::Z);
    vector<int> pyr_levels = pyr_info_levels();
    for (int i = 0; i < n_pyr_infos(); ++i)
    {
        int pyr_level = pyr_levels[i];
        if (pyr_level == 0)
            MCP3D_ASSERT(zratios.at(pyr_level) == 1)
        else
            MCP3D_ASSERT(zratios.at(pyr_level) >= 1)
        if (i > 0)
        {
            bool ratio_valid = false;
            int delta_pyr_levels = pyr_levels[i] - pyr_levels[i - 1];
            for (int j = 0; j <= delta_pyr_levels; ++j)
                if (pyr_levels[i - 1] * mcp3d::IntPow<int>(2, j) == pyr_levels[i])
                    ratio_valid = true;
            MCP3D_ASSERT(ratio_valid)
        }
        MCP3D_ASSERT(zratios.at(pyr_level) <= zdim(pyr_level))
    }
}

void mcp3d::MChannelInfo::RemovePyrInfoWithUndersizedVolumes()
{
    if (!HasPyrInfo(0))
        return;
    if (n_pyr_infos() == 1)
        return;
    for (const auto& axis: vector<mcp3d::ChannelAxis>({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X}))
    {
        vector<int> levels = pyr_info_levels();
        unordered_map<int, int> ratios = pyr_ratios(axis);
        for (const auto pyr_level: levels)
        {
            if (ratios.at(pyr_level) > channel_pyr_infos_.at(pyr_level).dim(axis))
            {
                RemovePyrInfo(pyr_level);
                MCP3D_MESSAGE("discarding MChannelPyrInfo at level " + to_string(pyr_level) +
                              " due to volume being under-sized along " + mcp3d::ChannelAxisStr(axis) + " axis")
            }
        }
    }
}
