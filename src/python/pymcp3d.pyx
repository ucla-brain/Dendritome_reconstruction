import os
from copy import deepcopy
import numpy as np
# np.PyArray_SimpleNewFromData gives seg fault without this line
np.import_array()
cimport numpy as np
from libc.stdint cimport uint8_t
from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from cython.operator cimport dereference as deref
from cython.operator cimport address as addr
from cpython.object cimport Py_EQ
from cython.view cimport array as cvarray
cimport mcp3d_cython

# VoxelType enum values
M8U = 0
M8S = 1
M16U = 2
M16S = 4
M32U = 5
M32S = 6
M64U = 7
M64S = 8
M32F = 9
M64F = 10
VOXEL_UNKNOWN = -1

# mapping from python int to mcp3d::VoxelType
cdef unordered_map[int, mcp3d_cython.VoxelType] _VoxelTypeIntEnumMapper
_VoxelTypeIntEnumMapper[M8U] = mcp3d_cython.M8U
_VoxelTypeIntEnumMapper[M8S] = mcp3d_cython.M8S
_VoxelTypeIntEnumMapper[M16U] = mcp3d_cython.M16U
_VoxelTypeIntEnumMapper[M16S] = mcp3d_cython.M16S
_VoxelTypeIntEnumMapper[M32U] = mcp3d_cython.M32U
_VoxelTypeIntEnumMapper[M32S] = mcp3d_cython.M32S
_VoxelTypeIntEnumMapper[M64U] = mcp3d_cython.M64U
_VoxelTypeIntEnumMapper[M64S] = mcp3d_cython.M64S
_VoxelTypeIntEnumMapper[M32F] = mcp3d_cython.M32F
_VoxelTypeIntEnumMapper[M64F] = mcp3d_cython.M64F
_VoxelTypeIntEnumMapper[VOXEL_UNKNOWN] = mcp3d_cython.VOXEL_UNKNOWN

# mapping from python int (type cast from VoxelType enum) to npy type
cdef unordered_map[int, int] _VoxelTypeNpyTypeEnumMapper
_VoxelTypeNpyTypeEnumMapper[M8U] = np.NPY_UINT8
_VoxelTypeNpyTypeEnumMapper[M8S] = np.NPY_INT8
_VoxelTypeNpyTypeEnumMapper[M16U] = np.NPY_UINT16
_VoxelTypeNpyTypeEnumMapper[M16S] = np.NPY_INT16
_VoxelTypeNpyTypeEnumMapper[M32U] = np.NPY_UINT32
_VoxelTypeNpyTypeEnumMapper[M32S] = np.NPY_INT32
_VoxelTypeNpyTypeEnumMapper[M64U] = np.NPY_UINT64
_VoxelTypeNpyTypeEnumMapper[M64S] = np.NPY_INT64
_VoxelTypeNpyTypeEnumMapper[M32F] = np.NPY_FLOAT32
_VoxelTypeNpyTypeEnumMapper[M64F] = np.NPY_FLOAT64

# FileFormat enum values
TIFF = 0
OMETIFF = 1
IMARIS = 2
HDF5 = 3
FILE_UNKNOWN = -1

# mapping from python int to mcp3d::FileFormat
cdef unordered_map[int, mcp3d_cython.FileFormat] _FileFormatsIntEnumMapper
_FileFormatsIntEnumMapper[TIFF] = mcp3d_cython.TIFF
_FileFormatsIntEnumMapper[OMETIFF] = mcp3d_cython.OMETIFF
_FileFormatsIntEnumMapper[IMARIS] = mcp3d_cython.IMARIS
_FileFormatsIntEnumMapper[HDF5] = mcp3d_cython.HDF5
_FileFormatsIntEnumMapper[FILE_UNKNOWN] = mcp3d_cython.FILE_UNKNOWN

# wrapper classes generally owns a pointer to the wrapped c++ instance.
# exposed python API calls c++ functions via the pointer to modify the
# wrapped instance. suppose python class A wraps c++ class B, which has as
# its member c++ class C instance. B provides getter to C. python class D
# wraps C. the corresponding getter in A requesting class D instance
# shall first instantiate default instance D, followed by copying D's
# pointed class C c++ object from B's c++ getter
# a = A()
# d = a.getC()
# type(d) = D
# class A:
#    cdef B* wrapped
#    ...
#    def getC(self):
#        d = D()
#        d.wrapped[0] = a.wrapped.getC()
#        return d
# class D:
#    cdef C* wrapped
#    ...
#
# A should not maintain states of B via its own python attributes.
# faithfully translate calls instead


cdef class MChannelInfo:
    cdef mcp3d_cython.MChannelInfo* wrapped

    # cdef class constructor essentially calls different overloads of
    # c++ constructor for wrapped class
    def __cinit__(self, image_info_path=None):
        if image_info_path is not None:
            self.wrapped = new mcp3d_cython.MChannelInfo(image_info_path)
        else:
            self.wrapped = new mcp3d_cython.MChannelInfo()

    def __iadd__(self, other):
        if not isinstance(other, MChannelInfo):
            raise TypeError('+= operator can only be used with MChannelInfo instances')
        self.wrapped.AddPyrInfos((<MChannelInfo>other).wrapped[0].channel_pyr_infos())

    def LevelPathsExist(self, int pyr_level):
        self.wrapped.LevelPathsExist(pyr_level)

    def ValidatePyramidsStructure(self):
        self.wrapped.Validate()

    def Save(self):
        self.wrapped.Save()

    def Clear(self):
        self.wrapped.Clear()

    def VoxelBytes(self, int pyr_level=0):
        return self.wrapped.VoxelBytes(pyr_level)

    def ImagePath(self, int pyr_level, int z, int y=0, int x=0):
        # self.wrapped.ImagePath(pyr_level, z, y, x)
        return self.wrapped.ImagePath(pyr_level, z, y, x).decode(encoding='ascii')

    def empty(self):
        return self.wrapped.empty()

    def n_pyr_levels(self):
        return self.wrapped.n_pyr_levels()

    def xyz_dims(self, int pyr_level=0):
        return self.wrapped.xyz_dims(pyr_level)

    def xdim(self, int pyr_level=0):
        return self.wrapped.xdim(pyr_level)

    def ydim(self, int pyr_level=0):
        return self.wrapped.ydim(pyr_level)

    def zdim(self, int pyr_level=0):
        return self.wrapped.zdim(pyr_level)

    def voxel_type(self, int pyr_level=0):
        return <int>self.wrapped.voxel_type(pyr_level)

    def file_format(self, int pyr_level=0):
        return <int>self.wrapped.file_format(pyr_level)

    def volume_root_dir(self):
        return self.wrapped.channel_root_dir()

    def pyr_xy_ratios(self):
        return self.wrapped.pyr_xy_ratios()

    def pyr_z_ratios(self):
        return self.wrapped.pyr_z_ratios()

    def __dealloc__(self):
        del self.wrapped


cdef class MImageInfo:
    cdef mcp3d_cython.MImageInfo* wrapped

    def __cinit__(self, other=None):
        self.wrapped = new mcp3d_cython.MImageInfo()

    def channel_info(self, int channel_number, int time=0):
        channel_info = MChannelInfo()
        channel_info.wrapped[0] = self.wrapped.channel_info(channel_number, time)
        return channel_info


cdef class MImageBlock:
    cdef mcp3d_cython.MImageBlock* wrapped

    def __cinit__(self, offsets=None, extents=None, strides=None):
        if not offsets:
            offsets = []
        if not extents:
            extents = []
        if not strides:
            strides = []
        self.wrapped = new mcp3d_cython.MImageBlock(offsets, extents, strides)

    def Clear(self):
        self.wrapped.Clear()

    def empty(self):
        return self.wrapped.empty()

    def offsets(self):
        return self.wrapped.offsets()

    def extents(self):
        return self.wrapped.extents()

    def strides(self):
        return self.wrapped.strides()

    def __richcmp__(self, other, Py_EQ):
        if not isinstance(other, MImageBlock):
            raise TypeError('equality can only be tested on two MImageBlock instances')
        return self.offsets() == other.offsets() and \
               self.extents() == other.extents() and \
               self.strides() == other.strides()

    def __dealloc__(self):
        del self.wrapped


cdef class MImageView:
    cdef mcp3d_cython.MImageView* wrapped

    def __cinit__(self):
        self.wrapped = new mcp3d_cython.MImageView()

    def Clear(self):
        self.wrapped.Clear()

    def PrintView(self):
        self.wrapped.PrintView()

    def ViewMemorySize(self, string unit):
        return self.wrapped.ViewMemorySize(unit)

    def pyr_level(self):
        return self.wrapped.pyr_level()

    def pyr_level_offsets(self, int channel_number, int time=0):
        return self.wrapped.pyr_level_offsets(channel_number, time)

    def pyr_level_extents(self, int channel_number, int time=0):
        return self.wrapped.pyr_level_extents(channel_number, time)

    def pyr_level_strides(self, int channel_number, int time=0):
        return self.wrapped.pyr_level_strides(channel_number, time)

    def view_channels(self):
        return self.wrapped.view_channels()

    def view_times(self):
        return self.wrapped.view_times()

    def empty(self):
        return self.wrapped.empty()

    def is_unit_strided(self, int channel_number, const string& axes):
        return self.wrapped.is_unit_strided(channel_number, axes)

    def dims(self):
        return self.wrapped.dims()

    def xyz_dims(self, int channel_number, int time=0):
        return self.wrapped.xyz_dims(channel_number, time)

    def xdim(self, int channel_number, int time=0):
        return self.wrapped.xdim(channel_number, time)

    def ydim(self, int channel_number, int time=0):
        return self.wrapped.ydim(channel_number, time)

    def zdim(self, int channel_number, int time=0):
        return self.wrapped.zdim(channel_number, time)

    def n_channels(self):
        return self.wrapped.n_channels()

    def n_times(self):
        return self.wrapped.n_times()

    def n_volumes(self):
        return self.wrapped.n_volumes()

    def VolumeVoxels(self, int channel_number, int time=0):
        return self.wrapped.VolumeVoxels(channel_number, time)

    def PlaneVoxels(self, int channel_number, int time=0):
        return self.wrapped.PlaneVoxels(channel_number, time)

    def VoxelBytes(self):
        return self.wrapped.VoxelBytes()

    def VolumeBytes(self, int channel_number, int time=0):
        return self.wrapped.VolumeBytes(channel_number, time)

    def PlaneBytes(self, int channel_number, int time=0):
        return self.wrapped.PlaneBytes(channel_number, time)

    def voxel_type(self):
        return <int>self.wrapped.voxel_type()

    def __dealloc__(self):
        del self.wrapped


cdef class MImageBase:
    # not available in Python-space:
    cdef mcp3d_cython.MImageBase* base

    # available in Python-space:
    cdef public object _data, _wrapped_data, _older_data

    def __cinit__(self):
        # do not new MImageBase() here, MImageBase is pure virtual c++ class
        # self.base, copied from self.wrapped in derived classes __cinit__
        # will be used in calls
        self._data = []
        self._wrapped_data = []
        self._older_data = []

    # return a MImageInfo instance with c++ MImageInfo instance copy constructed
    # from self.base.image_info().
    def image_info(self):
        info = MImageInfo()
        # copy contruct c++ MImageInfo instance wrapped in pointer
        info.wrapped[0] = self.base.image_info()
        return info

    def selected_view(self):
        view = MImageView()
        view.wrapped[0] = self.base.selected_view()
        return view

    def loaded_view(self):
        view = MImageView()
        view.wrapped[0] = self.base.loaded_view()
        return view

    def xyz_dims(self, int channel=0, int pyr_level=0, int time=0):
        return self.base.xyz_dims(channel, pyr_level, time)

    def dims(self):
        return self.base.dims()

    def n_volumes(self):
        return self.base.n_volumes()

    def loaded_pyr_level(self):
        return self.base.loaded_pyr_level()

    def selected_pyr_level(self):
        return self.base.selected_pyr_level()

    def n_pyr_levels(self, int channel=0, int time=0):
        return self.base.n_pyr_infos(channel, time)

    def __dealloc__(self):
         pass


cdef class MImage(MImageBase):
    cdef mcp3d_cython.MImage* wrapped

    def __cinit__(self, volume_root_path, channel_dir_names=None):
        if not os.path.exists(volume_root_path):
            raise ValueError('{} does not exist'.format(volume_root_path))
        # gaurantee byte array type for c++ string
        if isinstance(volume_root_path, str):
            volume_root_path = volume_root_path.encode()
        if channel_dir_names is None:
            channel_dir_names = []
        channel_dir_names = [channel_dir_name.encode()
                             if isinstance(channel_dir_name, str)
                             else channel_dir_name
                             for channel_dir_name in channel_dir_names]
        self.wrapped = new mcp3d_cython.MImage(volume_root_path, <vector[string]>channel_dir_names)
        self.base = self.wrapped

    def ReadImageInfo(self, channel_numbers=None, bool ignore_saved=False):
        # if channel numbers is None, assume reading channel 0
        if channel_numbers is None:
            channel_numbers = [0]
        elif isinstance(channel_numbers, int):
            channel_numbers = [channel_numbers]
        else:
            for channel_number in channel_numbers:
                if not isinstance(channel_number, int) or channel_number < 0:
                    raise ValueError('channel numbers must be non negative integer')
        # call wrapped instance method
        self.wrapped.ReadImageInfo(<vector[int]>channel_numbers, ignore_saved)

    def SaveImageInfo(self):
        self.wrapped.SaveImageInfo()

    def SelectView(self, MImageBlock block, channel_numbers, rl=0,
                   bool interpret_view_as_local=False):
        if isinstance(channel_numbers, int):
            channel_numbers = [channel_numbers]
        elif not isinstance(channel_numbers, list):
            raise ValueError('channel_numbers should be a list of integer or a scalar integer')
        self.wrapped.SelectView(block.wrapped[0], <vector[int]>channel_numbers,
                                rl, interpret_view_as_local)

    # combine ReadData and Volume() function call of c++ to read and return
    # selected data volume
    def ReadData(self, black_background=True, mode="verbose"):
        # c++ wrapped instance ReadData() call
        if isinstance(mode, str):
            mode = mode.encode()
        self._ReadData(black_background=black_background, mode=mode)

        # create numpy array wrapper
        cdef uint8_t* ptr
        cdef np.npy_intp view_dims[3]

        if self.wrapped.loaded_view().empty():
            return None
        ndims = 3
        # numpy axis order correspond to [z, y, x], with the fastest altering axis appearing last
        view_dims[0] = self.wrapped.loaded_view().zdim(0, 0)
        view_dims[1] = self.wrapped.loaded_view().ydim(0, 0)
        view_dims[2] = self.wrapped.loaded_view().xdim(0, 0)
        # wrap data in c++ class with numpy array
        data = [None] * self.wrapped.loaded_view().n_channels()
        for i in range(self.wrapped.loaded_view().n_channels()):
            ptr = self.base.Volume[uint8_t](i)
            npy_type = _VoxelTypeNpyTypeEnumMapper[<int>(self.wrapped.voxel_type())]
            # not sure how to avoid the copy here. setting NPY_OWNDATA does not work
            data[i] = np.PyArray_SimpleNewFromData(ndims, view_dims, npy_type, ptr).copy()
        # release data ownership from c++ object
        self.wrapped.ClearLoaded()
        if len(data) == 1:
            return data[0]
        return data

    def _ReadData(self, bool black_background=True, const string& mode=b"verbose"):
        self.wrapped.ReadData(black_background, mode)

    def WriteImagePyramids(self, int channel_number, int start_parent_level=0,
                           int end_parent_level=999, bool multi_threading=False,
                           bool save_image_info = True, int write_format=FILE_UNKNOWN):
        cdef mcp3d_cython.FileFormat write_format_c

        if _FileFormatsIntEnumMapper.find(write_format) == _FileFormatsIntEnumMapper.end():
            raise ValueError('write_format = {} is undefined'.format(write_format))
        write_format_c = _FileFormatsIntEnumMapper[write_format]
        self.wrapped.WriteImagePyramids(channel_number, start_parent_level,
                                        end_parent_level, multi_threading,
                                        save_image_info, write_format_c)

    def __dealloc__(self):
        del self.wrapped


