from libc.stdint cimport uint8_t
from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "image/mcp3d_voxel_types.hpp" namespace "mcp3d":
    cdef cppclass VoxelType:
        pass
    cdef VoxelType M8U "mcp3d::VoxelType::M8U"
    cdef VoxelType M8S "mcp3d::VoxelType::M8S"
    cdef VoxelType M16U "mcp3d::VoxelType::M16U"
    cdef VoxelType M16S "mcp3d::VoxelType::M16S"
    cdef VoxelType M32U "mcp3d::VoxelType::M32U"
    cdef VoxelType M32S "mcp3d::VoxelType::M32S"
    cdef VoxelType M64U "mcp3d::VoxelType::M64U"
    cdef VoxelType M64S "mcp3d::VoxelType::M64S"
    cdef VoxelType M32F "mcp3d::VoxelType::M32F"
    cdef VoxelType M64F "mcp3d::VoxelType::M64F"
    cdef VoxelType VOXEL_UNKNOWN "mcp3d::VoxelType::UNKNOWN"


cdef extern from "image/mcp3d_file_formats.hpp" namespace "mcp3d":
    cdef cppclass FileFormat:
        pass
    cdef FileFormat TIFF "mcp3d::FileFormat::TIFF"
    cdef FileFormat OMETIFF "mcp3d::FileFormat::OMETIFF"
    cdef FileFormat IMARIS "mcp3d::FileFormat::IMARIS"
    cdef FileFormat HDF5 "mcp3d::FileFormat::HDF5"
    cdef FileFormat FILE_UNKNOWN "mcp3d::FileFormat::UNKNOWN"


cdef extern from "image/mcp3d_channel_info.hpp" namespace "mcp3d":
    cdef cppclass MChannelPyrInfo:
        MChannelPyrInfo() except +
        MChannelPyrInfo(const string& json_path, FileFormat expected_format) except +
        MChannelPyrInfo(const MChannelPyrInfo& other);
        bool_t is_level0() except +
        int rank() except +
        int xdim() except +
        int ydim() except +
        int zdim() except +
        int n_channels() except +
        int n_times() except +
        FileFormat file_format() except +
        VoxelType voxel_type() except +
        vector[int] xyz_dims() except +
        vector[int] dims() except +
        int chunk_xdim() except +
        int chunk_ydim() except +
        int chunk_zdim() except +
        vector[int] chunk_xyz_dims() except +
        int n_xchunks() except +
        int n_ychunks() except +
        int n_zchunks() except +
        int n_total_chunks()  except +
        string channel_pyr_dir() except +
        const string &dims_order() except +
        const vector[string]& slice_image_names() except +


# do not supply default argument. it causes the compiler to complain about
# incomplete type
cdef extern from "image/mcp3d_channel_info.hpp" namespace "mcp3d":
    cdef cppclass MChannelInfo:
        MChannelInfo() except +
        MChannelInfo(const string& image_info_path) except +
        void AddPyrInfos(vector[MChannelPyrInfo]& channel_pyr_infos) except +
        void AddPyrInfo(MChannelPyrInfo& pyr_info) except +
        bool_t LevelPathsExist(int pyr_level) except +
        void Validate() except +
        void Save() except +
        void Clear() except +
        string ImagePath(int pyr_level, int z, int y, int x) except +
        bool_t empty() except +
        int n_pyr_levels() except +
        int VoxelBytes(int pyr_level) except +
        vector[int] xyz_dims(int pyr_level) except +
        int xdim(int pyr_level) except +
        int ydim(int pyr_level) except +
        int zdim(int pyr_level) except +
        VoxelType voxel_type(int pyr_level) except +
        FileFormat file_format(int pyr_level) except +
        string channel_root_dir() except +
        vector[MChannelPyrInfo]& channel_pyr_infos() except +
        const vector[int] pyr_xy_ratios() except +
        const vector[int] pyr_z_ratios() except +


cdef extern from "image/mcp3d_image_info.hpp" namespace "mcp3d":
    cdef cppclass MImageInfo:
        MImageInfo() except +
        MImageInfo(const string& volume_root_dir,
                   const vector[string]& channel_dir_names) except +
        void Save() except +
        void Clear() except +
        const MChannelInfo& channel_info(int channel, int time) except +


cdef extern from "image/mcp3d_image_view.hpp" namespace "mcp3d":
    cdef cppclass MImageBlock:
        MImageBlock(const vector[int] offsets, const vector[int] extents,
                    const vector[int] strides) except +
        bool_t operator== (const MImageBlock& other) except +
        bool_t empty() except +
        void Clear() except +
        vector[int] offsets() except +
        vector[int] extents() except +
        vector[int] strides() except +


cdef extern from "image/mcp3d_image_view.hpp" namespace "mcp3d":
    cdef cppclass MImageView:
        MImageView() except +
        MImageView(const MChannelInfo& img_info) except +
        void SelectView(const MImageBlock &image_block, int pyr_level,
                        bool_t interpret_view_as_local, VoxelType voxel_type) except +
        void Clear() except +
        void PrintView() except +
        double ViewMemorySize(const string &unit) except +
        int pyr_level()  except +
        const vector[int]& pyr_level_offsets(int channel_number, int time) except +
        const vector[int]& pyr_level_extents(int channel_number, int time) except +
        const vector[int]& pyr_level_strides(int channel_number, int time) except +
        const vector[int]& view_channels() except +
        const vector[int]& view_times() except +
        bool_t empty() except +
        bool_t is_unit_strided(int channel_number, const string& axes) except +
        vector[int] xyz_dims(int channel, int time) except +
        vector[vector[int]] dims() except +
        int xdim(int channel_number, int time) except +
        int ydim(int channel_number, int time) except +
        int zdim(int channel_number, int time) except +
        int n_channels() except +
        int n_times() except +
        int n_volumes() except +
        long VolumeVoxels(int channel_number, int time) except +
        long PlaneVoxels(int channel_number, int time) except +
        int VoxelBytes() except +
        long VolumeBytes(int channel_number, int time) except +
        long PlaneBytes(int channel_number, int time) except +
        VoxelType voxel_type() except +


cdef extern from "image/mcp3d_image_base.hpp" namespace "mcp3d":
    cdef cppclass MImageBase:
        MImageBase() except +
        MImageBase(bool_t from_disk) except +
        uint8_t* Volume[uint8_t](int c, int t) except +
        uint8_t* Volume[uint8_t](int index) except +
        uint8_t* ConstVolume[uint8_t] (int volume_index) except +
        const MImageInfo& image_info() except +
        const MImageView& selected_view() except +
        const MImageView& loaded_view() except +
        const vector[MChannelPyrInfo]& image_pyr_infos() except +
        vector[int] xyz_dims(int channel, int pyr_level, int time) except +
        vector[int] dims() except +
        int selected_pyr_level() except+
        int loaded_pyr_level() except+
        int n_pyr_infos(int channel, int time) except +
        int n_volumes() except +
        VoxelType voxel_type() except +


cdef extern from "image/mcp3d_image.hpp" namespace "mcp3d":
    cdef cppclass MImage(MImageBase):
        MImage() except +
        MImage(const string &img_root_dir, const vector[string]& channel_dir_names) except +
        void ReadImageInfo(const vector[int]& channel_numbers, bool_t ignore_saved) except +
        void ReadImageInfo(int channel_number, bool_t ignore_saved) except +
        void SaveImageInfo() except +
        void SelectView(const MImageBlock &global_view, const vector[int]& channel_numbers,
                        int rl, bool_t interpret_view_as_local) except +
        void ReadData(bool_t black_background, const string& mode) except +
        void ClearLoaded() except +
        void WriteImagePyramids(int channel_number, int start_parent_level,
                                int end_parent_level, bool_t multi_threading,
                                bool_t save_image_info, FileFormat write_format) except +
        void ReleaseData() except +


