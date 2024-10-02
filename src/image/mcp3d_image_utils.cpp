//
// Created by muyezhu on 10/13/17.
//
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <boost/filesystem.hpp>
#include <opencv2/core/core.hpp>
#include "common/mcp3d_utility.hpp"
#include "mcp3d_image_utils.hpp"

using namespace std;

int64_t mcp3d::LinearAddressRobust(const std::vector<int> &dims_, int d0, int d1, int d2)
{
    MCP3D_ASSERT(dims_.size() >= 2 && dims_.size() <= 3)
    vector<int> dims(dims_);
    if (dims.size() == 2)
        dims.push_back(1);
    MCP3D_ASSERT(d0 >= 0 && d0 < dims[0] &&
                 d1 >= 0 && d1 < dims[1] &&
                 d2 >= 0 && d2 < dims[2])
    std::vector<int64_t> strides(3, 0);
    for (int i = 2; i >= 0; --i)
    {
        // if dims[i] == 0, dimension stride is 0
        if (dims[i] > 1)
        {
            int64_t s = 1;
            for (int j = 2; j > i; --j)
                s *= (int64_t)dims[j];
            strides[i] = s;
        }
    }
    return (int64_t)d0 * strides[0] + (int64_t)d1 * strides[1] + (int64_t)d2 * strides[2];
}

int64_t mcp3d::LinearAddress(const std::vector<int> &dims_, int d0, int d1, int d2)
{
    MCP3D_ASSERT(dims_.size() >= 2 && dims_.size() <= 3)
    vector<int> dims(dims_);
    if (dims.size() == 2)
        dims.push_back(1);
    std::vector<int64_t> strides(3, 0);
    for (int i = 2; i >= 0; --i)
    {
        if (dims[i] > 1)
        {
            int64_t s = 1;
            for (int j = 2; j > i; --j)
                s *= (int64_t)dims[j];
            strides[i] = s;
        }
    }
    return (int64_t)d0 * strides[0] + (int64_t)d1 * strides[1] + (int64_t)d2 * strides[2];
}

array<int, 3> mcp3d::ZyxFromLinearAddress(const vector<int> &dims, int64_t address)
{
    int zyx_ptr[3];
    mcp3d::ZyxFromLinearAddress(dims, address, zyx_ptr);
    array<int, 3> zyx = {zyx_ptr[0], zyx_ptr[1], zyx_ptr[2]};
    return zyx;
}

void mcp3d::ZyxFromLinearAddress(const std::vector<int> &dims, int64_t address, int *zyx)
{
    MCP3D_ASSERT(dims.size() == (size_t)3);
    MCP3D_ASSERT(dims[0] > 0 && dims[1] > 0 && dims[2] > 0);
    if (address < 0 || address >= mcp3d::ReduceProdSeq<int64_t>(dims))
        MCP3D_OUT_OF_RANGE("address out of range for volume with dimensions " + mcp3d::JoinVector(dims))
    MCP3D_ASSERT(zyx)
    int ydim = dims[1], xdim = dims[2];
    zyx[0] = (int)(address / ((int64_t)ydim * (int64_t)xdim));
    zyx[1] = (int)(address % ((int64_t)ydim * (int64_t)xdim) / (int64_t)xdim);
    zyx[2] = (int)(address % (int64_t)xdim);
}

void mcp3d::NeighborAddresses(const vector<int> &dims, int64_t address, int connectivity_type,
                              unique_ptr<int64_t []>& neighbor_addresses)
{
    int zyx[3];
    mcp3d::ZyxFromLinearAddress(dims, address, zyx);
    int z = zyx[0], y = zyx[1], x = zyx[2];
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
}

void mcp3d::NeighborAddresses(const vector<int> &dims, int z, int y, int x,
                              int connectivity_type, unique_ptr<int64_t[]>& neighbor_addresses)
{
    MCP3D_ASSERT(dims.size() == (size_t)3)
    MCP3D_ASSERT(connectivity_type == 1 || connectivity_type == 2 || connectivity_type == 3)
    int64_t address = mcp3d::LinearAddressRobust(dims, z, y, x);
    neighbor_addresses = make_unique<int64_t []>(27);
    auto neighbor_addresses_ptr = neighbor_addresses.get();
    int64_t dims_[3] = {(int64_t)dims[0], (int64_t)dims[1], (int64_t)dims[2]};
    int64_t n_neighbors = 0, n_yx_voxels = (int64_t)dims[1] * (int64_t)dims[2],
            z_ = (int64_t)z, y_ = (int64_t)y, x_ = (int64_t)x;
    int64_t z_start = std::max(0l, z_ - 1), y_start = std::max(0l, y_ - 1), x_start = std::max(0l, x_ - 1),
            z_end = std::min(dims_[0], z_ + 2), y_end = std::min(dims_[1], y_ + 2), x_end = std::min(dims_[2], x_ + 2);
    int64_t dz_abs, dy_abs, z_offset, y_offset;
    for (int64_t nz = z_start; nz < z_end; ++nz)
    {
        dz_abs = abs(nz - z_);
        z_offset = nz * n_yx_voxels;
        for (int64_t ny = y_start; ny < y_end; ++ny)
        {
            dy_abs = abs(ny - y_);
            y_offset = ny * dims_[1];
            for (int64_t nx = x_start; nx < x_end; ++nx)
                // filter connectivity type
                if (dz_abs + dy_abs + abs(nx - x_) <= connectivity_type)
                {
                    int64_t neighbor_address = z_offset + y_offset + nx;
                    if (neighbor_address != address)
                        neighbor_addresses_ptr[++n_neighbors] = neighbor_address;
                }
        }
    }
    neighbor_addresses_ptr[0] = n_neighbors;

}

int mcp3d::ResidualStride(int extent, int stride)
{
    MCP3D_ASSERT(extent > 0 && stride > 0)
    int n_in_boundary_voxels = mcp3d::StridedExtent(extent, stride);
    int unit_stride_extent = n_in_boundary_voxels * stride;
    return unit_stride_extent <= extent ? 0 : unit_stride_extent - extent;
}

bool mcp3d::VolumeEqual(const void *ptr1, const void *ptr2, const vector<int> &dims, int element_size)
{
    MCP3D_ASSERT(element_size > 0)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    size_t n_bytes = mcp3d::ReduceProdSeq<size_t>(dims) * element_size;
    return memcmp(ptr1, ptr2, n_bytes) == 0;
}

void mcp3d::CopyDataVolume(void *dst, const vector<int> &dst_dims,
                           void *src, const vector<int> &src_dims, int type_size,
                           const MImageBlock &dst_block, const MImageBlock &src_block)
{
    MCP3D_ASSERT(dst && src)
    MCP3D_ASSERT(dst_dims.size() == 3 && src_dims.size() == 3)
    MCP3D_ASSERT(AllPositive(dst_dims) && AllPositive(src_dims))
    MCP3D_ASSERT(type_size > 0)
    std::vector<int> src_offsets(move(src_block.offsets())),
                     src_extents(move(src_block.extents())),
                     src_strides(move(src_block.strides())),
                     dst_offsets(move(dst_block.offsets())),
                     dst_strides(move(dst_block.strides()));
    // for any missing extents value (0 value), fill it in as if
    // offsets: end of dimension was given
    for (int i = 0; i < 3; ++i)
        if (src_extents[i] == 0)
            src_extents[i] = src_dims[i] - src_offsets[i];

    MCP3D_ASSERT(AllNonNegative(src_offsets) && AllNonNegative(src_extents))
    MCP3D_ASSERT(AllPositive(src_strides))
    MCP3D_ASSERT(AllNonNegative(SubtractSeq<int>(src_dims, AddSeq<int>(src_offsets, src_extents))))
    if (ReduceProdSeq<int>(src_extents) == 0)
        return;

    MCP3D_ASSERT(AllNonNegative(dst_offsets))
    MCP3D_ASSERT(AllPositive(dst_strides))
    std::vector<int> src_extents_strided =
            mcp3d::StridedExtents(src_extents, src_strides);
    std::vector<int> dst_extents = mcp3d::MultiplySeq<int>(dst_strides,
                                                           src_extents_strided);
    std::vector<int> dst_remain_extents = SubtractSeq<int>(dst_dims, AddSeq<int>(dst_offsets, dst_extents));
    MCP3D_ASSERT(AllNonNegative(dst_remain_extents))

    // treat the array as uint8_t type
    auto dst_uint8 = (uint8_t*)dst;
    auto src_uint8 = (uint8_t*)src;

    // copy data
    for (int z_dst = dst_offsets[0], z_src = src_offsets[0];
         z_src < src_offsets[0] + src_extents[0];
         z_dst += dst_strides[0], z_src += src_strides[0])
    {
        // memcpy for simple case
        if (mcp3d::ReduceProdSeq<int>(src_strides) == 1 &&
            mcp3d::ReduceProdSeq<int>(dst_strides) == 1)
        {
            if (mcp3d::AllZeros(src_offsets) && mcp3d::AllZeros(dst_offsets) &&
                src_extents == src_dims && dst_extents == dst_dims && src_dims == dst_dims)
            {
                #ifdef VERBOSE
                if (z_dst == dst_offsets[2])
                    MCP3D_MESSAGE("copy continuous block to same type")
                #endif
                size_t n_bytes = (size_t)type_size * mcp3d::ReduceProdSeq<int>(src_dims);
                memcpy(dst, src, n_bytes);
            }
            else
            {
                #ifdef VERBOSE
                if (z_dst == dst_offsets[2])
                    MCP3D_MESSAGE("copy discontinuous unit stride block to sampe type");
                #endif
                // element address. multiply by type size for byte address
                long addr_dst = mcp3d::LinearAddress(dst_dims, z_dst, dst_offsets[1], dst_offsets[2]),
                     addr_src = mcp3d::LinearAddress(src_dims, z_src, src_offsets[1], src_offsets[2]);
                for(int y_src = src_offsets[1]; y_src < src_offsets[1] + src_extents[1]; ++y_src)
                {
                    size_t n_bytes = src_extents[2] * (size_t)type_size;
                    memcpy(dst_uint8 + addr_dst * type_size,
                           src_uint8 + addr_src * type_size, n_bytes);
                    addr_dst += dst_dims[2];
                    addr_src += src_dims[2];
                }
            }
        }
        else
        {
            #ifdef VERBOSE
            if (z_dst == dst_offsets[2])
                MCP3D_MESSAGE("general copy")
            #endif
            for (int y_dst = dst_offsets[1], y_src = src_offsets[1];
                 y_src < src_offsets[1] + src_extents[1];
                 y_dst += dst_strides[1], y_src += src_strides[1])
            {
                int x_dst = dst_offsets[2], x_src = src_offsets[2];
                long addr_src = mcp3d::LinearAddress(src_dims, z_src, y_src, x_src),
                     addr_dst = mcp3d::LinearAddress(dst_dims, z_dst, y_dst, x_dst);
                for (; x_src < src_offsets[2] + src_extents[2]; x_src += src_strides[2])
                {
                    // copy type_size bytes at each element
                    for (int k = 0; k < type_size; ++k)
                        dst_uint8[addr_dst * type_size + k] = src_uint8[addr_src * type_size + k];
                    addr_dst += dst_strides[2];
                    addr_src += src_strides[2];
                }
            }
        }
    }
}

