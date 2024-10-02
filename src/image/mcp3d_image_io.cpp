//
// Created by muyezhu on 2/12/18.
//
#include <memory>
#include <algorithm>
#include <unordered_set>
#include "common/mcp3d_utility.hpp"
#include "mcp3d_image_macros.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_tiff_format.hpp"
#include "mcp3d_ome_tiff_format.hpp"
#include "mcp3d_imaris_format.hpp"
#include "mcp3d_hdf5_format.hpp"
#include "mcp3d_image_io.hpp"

using namespace std;

mcp3d::MImageIO::MImageIO()
{
    io_formats_[FileFormat::TIFF] = make_unique<mcp3d::MTiffFormat>();
    io_formats_[FileFormat::OMETIFF] = make_unique<mcp3d::MOmeTiffFormat>();
    io_formats_[FileFormat::IMARIS] = make_unique<mcp3d::MImarisFormat>();
}

void mcp3d::MImageIO::ReadChannelInfo(mcp3d::MImage& image, const string& channel_name, int pyr_level, bool load_from_json)
{
    try
    {
        if (!image.image_info().volume_layout().HasChannel(channel_name))
            MCP3D_RUNTIME_ERROR("channel name " + channel_name + " does not exist under volume " + image.image_info().volume_layout().volume_root_dir())
        cout << "Reading channel info from channel root directory " << image.image_info().volume_layout().channel_root_dir(channel_name, 0) << endl;
        image.image_info().channel_info(channel_name).RefreshChannelLayout();
        if (pyr_level > 0)
            MCP3D_ASSERT(image.image_info().channel_info(channel_name).HasPyrLevel(pyr_level))
        // if pyr_level < 0, all levels under channel_info.channel_layout() need to be read
        // else only level = pyr_level needs to be read
        vector<int> levels = pyr_level >= 0 ? vector<int>({pyr_level}) : image.image_info().channel_info(channel_name).channel_layout().pyr_levels();
        ReadChannelPyrInfos(image.image_info().channel_info(channel_name), levels, load_from_json);
        image.image_info().channel_info(channel_name).RemovePyrInfoWithUndersizedVolumes();
        image.image_info().channel_info(channel_name).AssertValid();
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

void mcp3d::MImageIO::ReadChannelPyrInfos(MChannelInfo& channel_info, vector<int> pyr_levels, bool load_from_json)
{
    try
    {
        for (const auto& pyr_level: pyr_levels)
            channel_info.RemovePyrInfo(pyr_level);
        // load existing __channel_info__.json
        if (load_from_json)
            channel_info.Load();
        for (auto& pyr_level: pyr_levels)
        {
            mcp3d::FileFormat pyr_level_format = channel_info.channel_layout().file_format(pyr_level);
            if (channel_info.HasPyrInfo(pyr_level))
            {
                if (channel_info.file_format(pyr_level) == pyr_level_format)
                    continue;
            }
            if (pyr_level_format == mcp3d::FileFormat::UNKNOWN)
            {
                MCP3D_MESSAGE("no known image format found under " + channel_info.channel_pyr_dir(pyr_level))
                continue;
            }
            MCP3D_ASSERT(io_formats_.find(pyr_level_format) != io_formats_.end())
            channel_info.AddPyrInfo(io_formats_.at(pyr_level_format)->ReadChannelPyrInfo(channel_info.channel_layout().channel_pyr_slices(pyr_level), pyr_level));
        }
    }
    catch(...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

// MChannelPyrInfo should have been read for selected channels at pyr_level
void mcp3d::MImageIO::ReadData(mcp3d::MImage &image)
{
    if (image.image_info().empty())
        MCP3D_RUNTIME_ERROR("MImageInfo is empty")
    if (image.selected_view().empty())
        MCP3D_RUNTIME_ERROR("no image view selected for reading")
    IMAGE_SELECTED_TYPED_CALL(mcp3d::MImageIO::FillViewWithZeros, image, image);
    if (image.selected_view().OutOfPyrImageBoundary())
    {
        #ifdef VERBOSE
        MCP3D_MESSAGE("view is entirely out of view level image bounds, fill with background voxels");
        #endif
        return;
    }
    int selected_level = image.selected_view().pyr_level();
    for (int view_channel: image.selected_view().view_channels())
    {
        MCP3D_ASSERT(image.image_info().HasChannelPyrInfo(view_channel, selected_level))
        mcp3d::FileFormat format = image.image_info().channel_info(view_channel, 0).file_format(selected_level);
        if (!ReadableFormat(format))
            MCP3D_DOMAIN_ERROR("reading " + mcp3d::FileFormatToString(format) + " format not supported")
        MCP3D_TRY(io_formats_[format]->ReadChannelData(image, view_channel);)
    }
}

void mcp3d::MImageIO::WriteViewVolume(const MImage &image, const string &out_dir, const string &img_name_prefix, FileFormat write_format)
{
    MCP3D_ASSERT(!image.loaded_view().empty())
    MCP3D_ASSERT(WritableFormat(write_format))
    mcp3d::MakeDirectories(out_dir);
    io_formats_[write_format]->WriteViewVolume(image, out_dir, img_name_prefix);
}

void mcp3d::MImageIO::WriteImagePyramid(mcp3d::MImage &image, int channel_number, int parent_level, bool multi_threading, mcp3d::FileFormat write_format)
{
    string pyr_level_dir;
    try
    {
        MCP3D_ASSERT(image.image_info().HasChannelPyrInfo(channel_number, parent_level))
        pyr_level_dir = mcp3d::JoinPath(image.image_info().volume_layout().channel_root_dir(channel_number, 0),
                                        mcp3d::MChannelLayout::PyrLevelDirName(parent_level + 1));
        MCP3D_ASSERT(mcp3d::RemovePath(pyr_level_dir))
        mcp3d::MakeDirectories(pyr_level_dir);
        mcp3d::FileFormat read_format = image.image_info().channel_info(channel_number, 0).file_format(parent_level);
        if (write_format == mcp3d::FileFormat::UNKNOWN)
            write_format = read_format;
        MCP3D_ASSERT(ReadableFormat(read_format))
        MCP3D_ASSERT(WritableFormat(write_format))
        io_formats_[write_format]->WriteImagePyramid(image, channel_number,
                                                     parent_level,
                                                     multi_threading);
    }
    catch (...)
    {
        MCP3D_MESSAGE("error during image pyramid creation. removing output directory " + pyr_level_dir)
        if (!mcp3d::RemovePath(pyr_level_dir))
            MCP3D_MESSAGE("cannot remove output directory " + pyr_level_dir)
        MCP3D_RETHROW(current_exception());
    }
}

bool mcp3d::MImageIO::ReadableFormat(FileFormat format)
{
    if (io_formats_.find(format) == io_formats_.end())
        return false;
    if (format == FileFormat::UNKNOWN)
        return false;
    return io_formats_[format]->CanRead();
}

bool mcp3d::MImageIO::WritableFormat(FileFormat format)
{
    if (io_formats_.find(format) == io_formats_.end())
        return false;
    if (format == FileFormat::UNKNOWN)
        return false;
    return io_formats_[format]->CanWrite();
}

#if MCP3D_MPI_BUILD

void mcp3d::MImageIO::WriteImagePyramidMPI(const std::string &img_root_dir, int channel,
                                           int parent_level, FileFormat write_format,
                                           bool abort_all_on_fail, MPI_Comm writer_comm)
{
    MCP3D_ASSERT(mcp3d::MPIInitialized())
    int rank, size;
    MPI_Comm_rank(writer_comm, &rank);
    MPI_Comm_size(writer_comm, &size);

    // set up logs. rank0 will only maintain log and make file system operations
    string log_dir(mcp3d::JoinPath({img_root_dir, "log"}));
    RANK0_COMM_CALL_SYNC(mcp3d::MakeDirectories, writer_comm, log_dir)
    time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    char time_cstr[13];
    strftime(time_cstr, 13, "%Y%m%d%H%M", localtime(&now));
    string err_log_path = mcp3d::JoinPath({log_dir, string("write_image_pyramid_err") + time_cstr}),
           out_log_path = mcp3d::JoinPath({log_dir, string("write_image_pyramid_out") + time_cstr});

    ofstream ofs;
    if (rank == 0)
    {
        ofs.open(err_log_path, ofstream::out);
        ofs.close();
        ofs.open(out_log_path, ofstream::out);
        ofs.close();
    }

    INIT_TIMER(0)
    TIC_COMM(writer_comm, 0)
    mcp3d::MImage image {};
    if (rank == 0)
    {
        image.ReadImageInfo(img_root_dir, channel);
        image.SaveImageInfo();
    }
    MPI_Barrier(writer_comm);
    if (rank != 0)
        image.ReadImageInfo(img_root_dir, channel);
    TOC_COMM(writer_comm, 0)

    if (rank == 0)
    {
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        ofs.open(out_log_path, ofstream::out | ofstream::app);
        ofs << "discarded " << world_size - size << " slow processes" << endl;
        ofs.close();

        if (image.image_info().empty())
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            ofs << "no known image data found under " << image.image_root_dir() << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (parent_level < 0 || parent_level >= image.n_pyr_levels())
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            ofs << "requested parent image level does not exist" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(writer_comm);
    mcp3d::FileFormat read_format = image.image_info().file_format(channel, parent_level);
    if (write_format == mcp3d::FileFormat::UNKNOWN)
        write_format = read_format;
    if (rank == 0)
    {
        ofs.open(out_log_path, ofstream::out | ofstream::app);
        ofs << "image root directory: " << img_root_dir << endl;
        ofs << "create one level of pyramid image from level " << parent_level << endl;
        ofs << "    input dimensions " << mcp3d::JoinVector<int>(image.xyz_dims(channel, parent_level), ", ") << endl;
        ofs << "    voxel type: " << mcp3d::VoxelTypeToString(
                image.voxel_type()) << endl;
        ofs << "reading image info: " << ELAPSE(0) << " seconds"<< endl;
        ofs.close();

        if(!ReadableFormat(read_format))
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            ofs << "image format unreadable" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (!WritableFormat(write_format))
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            ofs << "requested write format not supported" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(writer_comm);

    // remove older output, create directory
    string pyr_level_dir(mcp3d::JoinPath({image.channel_root_dir(channel, 0),
                                         mcp3d::MChannelLayout::PyrLevelDirName(parent_level + 1)}));
    RANK0_COMM_CALL_SYNC(mcp3d::RemovePath, writer_comm, pyr_level_dir)
    RANK0_COMM_CALL_SYNC(mcp3d::MakeDirectories, writer_comm, pyr_level_dir)

    string parent_dir(image.channel_pyr_info(channel, parent_level).channel_pyr_dir());
    int n_parent_imgs = image.channel_info(channel).zdim(parent_level);
    int n_workers = min(size - 1, n_parent_imgs);
    MPI_Status status;

    // check paths integrity
    if (rank > 0 && rank <= n_workers)
    {
        string img_path;
        int32_t img_id = 0;
        for (int i = (rank - 1) * n_parent_imgs / n_workers;
             i < min(rank * n_parent_imgs / n_workers, n_parent_imgs); ++i)
        {
            img_path = image.channel_info(channel).ImagePath(parent_level, i, 0,
                                                             0);
            img_id = i;
            // send failures
            if (!mcp3d::IsFile(img_path))
                MPI_Send(&img_id, 1, MPI_INT32_T, 0, 0, writer_comm);
        }
        // send -1 to mark task completion
        img_id = -1;
        MPI_Send(&img_id, 1, MPI_INT32_T, 0, 0, writer_comm);
    }
    if (rank == 0)
    {
        unordered_set<int32_t> bad_paths_id;
        int finished = 0;
        int32_t buffer;
        while (finished < n_workers)
        {
            MPI_Recv(&buffer, 1, MPI_INT32_T, MPI_ANY_SOURCE, MPI_ANY_TAG,
                     writer_comm, &status);
            if (buffer >= 0)
                bad_paths_id.insert(buffer);
            else
                ++finished;
        }
        if (!bad_paths_id.empty())
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            ofs << "bad image paths:" << endl;
            for (const auto& bad_path_id: bad_paths_id)
                ofs << image.channel_info(channel).ImagePath(parent_level,
                                                             bad_path_id, 0, 0) << endl;
            ofs <<"aborting" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    TIC_COMM(writer_comm, 0)
    io_formats_[write_format]->WriteImagePyramidMPI(image, channel, parent_level,
                                                    abort_all_on_fail,
                                                    err_log_path,
                                                    out_log_path, writer_comm);
    TOC_COMM(writer_comm, 0)
    if (rank == 0)
    {
        ofs.open(out_log_path, ofstream::out | ofstream::app);
        double n_min = ELAPSE(0) / 60;
        ofs << "image pyramid generation complete in " << n_min
            << " minutes with " << size - 1 << " writers" << endl;
        ofs.close();
    }
}

#endif



