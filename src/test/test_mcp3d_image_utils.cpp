//
// Created by muyezhu on 3/3/18.
//
#include <iostream>
#include <memory>
#include <cstring>
#include <algorithm>
#include <gtest/gtest.h>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_types.hpp"
#include "image/mcp3d_image_utils.hpp"

using namespace std;

TEST(ImageUtils, LinearAddressRobust)
{
    EXPECT_EQ(0, mcp3d::LinearAddressRobust({1, 4}, 0, 0));
    EXPECT_EQ(1, mcp3d::LinearAddressRobust({1, 4}, 0, 1));
    EXPECT_EQ(3, mcp3d::LinearAddressRobust({1, 4}, 0, 3));
    EXPECT_EQ(6, mcp3d::LinearAddressRobust({3, 4}, 1, 2));
    EXPECT_EQ(26, mcp3d::LinearAddress({3, 3, 3}, 2, 2, 2));
    EXPECT_EQ(1, mcp3d::LinearAddress({3, 1, 3}, 0, 0, 1));
    EXPECT_EQ(4, mcp3d::LinearAddress({3, 1, 3}, 1, 0, 1));
    EXPECT_EQ(4, mcp3d::LinearAddress({1, 3, 3}, 1, 1, 1));
    EXPECT_EQ(1, mcp3d::LinearAddress({1, 1, 3}, 1, 1, 1));
    EXPECT_EQ(234, mcp3d::LinearAddress({5, 10, 10}, 2, 3, 4));

    vector<vector<int>> dims_vec = {{1, 4}, {1}, {1, 2, 3, 4}};
    for (const auto& dims: dims_vec)
    {
        EXPECT_THROW(int64_t _ = mcp3d::LinearAddressRobust(dims, 1, 3),
                     mcp3d::MCPAssertionError);
    }
}

TEST(ImageUtils, LinearAddress)
{
    EXPECT_EQ(0, mcp3d::LinearAddress({1, 4}, 0, 0));
    EXPECT_EQ(1, mcp3d::LinearAddress({1, 4}, 0, 1));
    EXPECT_EQ(3, mcp3d::LinearAddress({1, 4}, 0, 3));
    EXPECT_EQ(6, mcp3d::LinearAddress({3, 4}, 1, 2));
    EXPECT_EQ(26, mcp3d::LinearAddress({3, 3, 3}, 2, 2, 2));
    EXPECT_EQ(1, mcp3d::LinearAddress({3, 1, 3}, 0, 0, 1));
    EXPECT_EQ(4, mcp3d::LinearAddress({3, 1, 3}, 1, 0, 1));
    EXPECT_EQ(4, mcp3d::LinearAddress({1, 3, 3}, 1, 1, 1));
    EXPECT_EQ(1, mcp3d::LinearAddress({1, 1, 3}, 1, 1, 1));
    EXPECT_EQ(234, mcp3d::LinearAddress({5, 10, 10}, 2, 3, 4));
}

TEST(ImageUtils, ZyxFromLinearAddress)
{
    vector<int> dims({256, 512, 1024});
    int z = 128, y = 45, x  = 912;
    int64_t address = mcp3d::LinearAddressRobust(dims, z, y, x);
    array<int, 3> zyx = mcp3d::ZyxFromLinearAddress(dims, address);
    EXPECT_EQ(z, zyx[0]);
    EXPECT_EQ(y, zyx[1]);
    EXPECT_EQ(x, zyx[2]);
    address = 0;
    zyx = mcp3d::ZyxFromLinearAddress(dims, address);
    EXPECT_EQ(0, zyx[0]);
    EXPECT_EQ(0, zyx[1]);
    EXPECT_EQ(0, zyx[2]);
    address = mcp3d::ReduceProdSeq<int64_t>(dims) - 1;
    zyx = mcp3d::ZyxFromLinearAddress(dims, address);
    EXPECT_EQ(dims[0] - 1, zyx[0]);
    EXPECT_EQ(dims[1] - 1, zyx[1]);
    EXPECT_EQ(dims[2] - 1, zyx[2]);
    address = -1;
    EXPECT_THROW(mcp3d::ZyxFromLinearAddress(dims, address), mcp3d::MCPOutOfRangeError);
    address = mcp3d::ReduceProdSeq<int64_t>(dims);
    EXPECT_THROW(mcp3d::ZyxFromLinearAddress(dims, address), mcp3d::MCPOutOfRangeError);
    EXPECT_THROW(mcp3d::ZyxFromLinearAddress({1, 2, 3, 4}, address), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::ZyxFromLinearAddress({1, 2, -3}, address), mcp3d::MCPAssertionError);
}

void ExpectedAddresses(const vector<int>& dims, int z, int y, int x,
                       int connectivity_type, unordered_set<int64_t>& expected_addresses)
{
    expected_addresses.clear();
    int64_t address = mcp3d::LinearAddress(dims, z, y, x);
    ASSERT_EQ((size_t)3, dims.size());
    ASSERT_TRUE(z >= 0 && z < dims[0]);
    ASSERT_TRUE(y >= 0 && z < dims[1]);
    ASSERT_TRUE(x >= 0 && z < dims[2]);

    int valid_zs[3] = {z, -1, -1};
    int valid_ys[3] = {y, -1, -1};
    int valid_xs[3] = {x, -1, -1};
    if (z - 1 >= 0)
        valid_zs[1] = z - 1;
    if (y - 1 >= 0)
        valid_ys[1] = y - 1;
    if (x - 1 >= 0)
        valid_xs[1] = x - 1;
    if (z + 1 < dims[0])
        valid_zs[2] = z + 1;
    if (y + 1 < dims[1])
        valid_ys[2] = y + 1;
    if (x + 1 < dims[2])
        valid_xs[2] = x + 1;
    int n_neighbors = 0;
    // if connectivity_type = 1 aka 6 neighbor, only one coordinates can be + or - 1
    for (int i = 1; i < 3; ++i)
    {
        if (valid_zs[i] >= 0)
            expected_addresses.insert(mcp3d::LinearAddress(dims, valid_zs[i], y, x));
        if (valid_ys[i] >= 0)
            expected_addresses.insert(mcp3d::LinearAddress(dims, z, valid_ys[i], x));
        if (valid_xs[i] >= 0)
            expected_addresses.insert(mcp3d::LinearAddress(dims, z, y, valid_xs[i]));
    }
    if (connectivity_type > 1)
    {
        // if connectivity_type = 2 aka 18 neighbor, two coordinates can be + or - 1
        for (int i = 1; i < 3; ++i)
            for (int j = 1; j < 3; ++j)
            {
                if (valid_ys[i] >= 0 && valid_xs[j] >= 0)
                    expected_addresses.insert(mcp3d::LinearAddress(dims, z, valid_ys[i], valid_xs[j]));
                if (valid_zs[i] >= 0 && valid_xs[j] >= 0)
                    expected_addresses.insert(mcp3d::LinearAddress(dims, valid_zs[i], y, valid_xs[j]));
                if (valid_zs[i] >= 0 && valid_ys[j] >= 0)
                    expected_addresses.insert(mcp3d::LinearAddress(dims, valid_zs[i], valid_ys[j], x));
            }
    }
    if (connectivity_type > 2)
    {
        // if connectivity_type = 3 aka 26 neighbor, all coordinates can be + or - 1
        for (int i = 1; i < 3; ++i)
            for (int j = 1; j < 3; ++j)
                for (int k = 1; k < 3; ++k)
                    if (valid_zs[i] >= 0 && valid_ys[j] >= 0 && valid_xs[k] >= 0)
                        expected_addresses.insert(mcp3d::LinearAddress(dims, valid_zs[i], valid_ys[j], valid_xs[k]));
    }
    expected_addresses.erase(address);
}

void VerifyNeighborAddresses(const unordered_set<int64_t>& expected_addresses,
                             const unique_ptr<int64_t []>& neighbor_addresses)
{
    EXPECT_EQ(expected_addresses.size(), (size_t)neighbor_addresses[0]);
    for (int i = 1; i < 27; ++i)
        if (i <= neighbor_addresses[0])
            EXPECT_TRUE(expected_addresses.find(neighbor_addresses[i]) != expected_addresses.end());
}

TEST(ImageUtils, NeighborAddresses)
{
    int zdim = 100, ydim = 100, xdim = 100, z, y, x;
    vector<int> dims({zdim, ydim, xdim});
    unique_ptr<int64_t []> neighbor_addresses;
    unordered_set<int64_t> expected_addresses;
    // connectivity 1
    int connectivity_type = 1;
    // no neighbor outside of range
    z = y = x = 50;
    INIT_TIMER(0)
    TIC(0)
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("NeighborAddresses", 0)
    EXPECT_EQ(6, neighbor_addresses[0]);
    TIC(0)
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("ExpectedAddresses", 0)
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // voxel on boundary along z axis
    z = 0;
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    EXPECT_EQ(5, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // connectivity 2
    connectivity_type = 2;
    // no neighbor outside of range
    z = 50;
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    EXPECT_EQ(18, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // voxel on boundary along y axis (5 neighbors eliminated)
    y = 99;
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    EXPECT_EQ(13, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // connectivity 3
    connectivity_type = 3;
    // no neighbor outside of range
    y = 50;
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    EXPECT_EQ(26, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // voxel on boundary along all axis. 7 neighbors remain
    z = y = x = 0;
    mcp3d::NeighborAddresses(dims, z, y, x, connectivity_type, neighbor_addresses);
    EXPECT_EQ(7, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
    // calling using address instead of z, y, x
    mcp3d::NeighborAddresses(dims, 0, connectivity_type, neighbor_addresses);
    EXPECT_EQ(7, neighbor_addresses[0]);
    ExpectedAddresses(dims, z, y, x, connectivity_type, expected_addresses);
    VerifyNeighborAddresses(expected_addresses, neighbor_addresses);
}

TEST(ImageUtils, ResidualStride)
{
    EXPECT_EQ(0, mcp3d::ResidualStride(4, 1));
    EXPECT_EQ(0, mcp3d::ResidualStride(4, 2));
    EXPECT_EQ(2, mcp3d::ResidualStride(4, 3));
    EXPECT_EQ(0, mcp3d::ResidualStride(4, 4));
    EXPECT_EQ(1, mcp3d::ResidualStride(4, 5));
    EXPECT_EQ(2, mcp3d::ResidualStride(4, 6));
    EXPECT_EQ(1, mcp3d::ResidualStride(5, 2));
}

TEST(ImageUtils, SetConstant)
{
    unique_ptr<uint16_t[]> data;
    vector<int> dims({10, 10, 10, 10});
    uint16_t c = 2515;
    EXPECT_THROW(mcp3d::SetConstant(data.get(), dims, c), mcp3d::MCPAssertionError);
    data = make_unique<uint16_t[]>(mcp3d::ReduceProdSeq<size_t>(dims));
    EXPECT_THROW(mcp3d::SetConstant(data.get(), {}, c), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::SetConstant(data.get(), {-1, 1}, c), mcp3d::MCPAssertionError);
    mcp3d::SetConstant(data.get(), dims, c);
    int addr_end = mcp3d::ReduceProdSeq<int>(dims);
    for (int i = 0; i < addr_end; ++i)
        EXPECT_EQ(c, data.get()[i]);
}

TEST(ImageUtils, SetConstantBlock)
{
    vector<int> dims({25, 25, 25}), dims_long ({25, 25, 25, 25}),
                offsets ({5, 4, 3}), extents ({9, 8, 7});
    unique_ptr<uint16_t[]> data_empty,
                           data = make_unique<uint16_t[]>(mcp3d::ReduceProdSeq<size_t>(dims));
    uint16_t* data_ptr = data.get();
    mcp3d::MImageBlock block_empty, block_oversize({0}, {30}),
                       block(offsets, extents);
    uint16_t background = 7, foreground = 77;
    long addr_end = mcp3d::ReduceProdSeq<long>(dims);
    // nullptr
    EXPECT_THROW(mcp3d::SetConstantBlock(data_empty.get(), dims, block_empty,
                                         background),
                 mcp3d::MCPAssertionError);
    // wrong dimension vector size
    EXPECT_THROW(mcp3d::SetConstantBlock(data.get(), dims_long, block_empty,
                                         background),
                 mcp3d::MCPAssertionError);
    // empty MImageBlock -> full black
    mcp3d::SetConstantBlock(data.get(), dims, block_empty, background);
    for (long i = 0; i < addr_end; ++i)
        EXPECT_EQ(background, data_ptr[i]);
    // oversize MImageBlock -> trunc to max value allowed along each dimension
    mcp3d::SetConstant(data.get(), dims, foreground);
    mcp3d::SetConstantBlock(data.get(), dims, block_oversize, background);
    for (long i = 0; i < addr_end; ++i)
        EXPECT_EQ(background, data_ptr[i]);
    // normal MImageBlock
    mcp3d::SetConstant(data.get(), dims, foreground);
    mcp3d::SetConstantBlock(data.get(), dims, block, background);
    for (int z = offsets[0]; z < offsets[0] + extents[0]; ++z)
        for (int y = offsets[1]; y < offsets[1] + extents[1]; ++y)
        {
            int x = offsets[2];
            int64_t addr = mcp3d::LinearAddressRobust(dims, z, y, x);
            for (; x < offsets[2] + extents[2]; ++x)
                EXPECT_EQ(background, data_ptr[addr++]);
        }
}

TEST(ImageUtils, CopyDataVolume)
{
    vector<int> dims_dst({5, 50, 50}), dims_dst_same({10, 100, 100}),
                dims_src ({10, 100, 100}),
                strides_fit({2, 2, 2}), strides_unit({1, 1, 1});
    unique_ptr<uint16_t[]> dst = make_unique<uint16_t[]>(mcp3d::ReduceProdSeq<size_t>(dims_dst)),
                           dst_same = make_unique<uint16_t[]>(mcp3d::ReduceProdSeq<size_t>(dims_dst_same)),
                           src = make_unique<uint16_t[]>(mcp3d::ReduceProdSeq<size_t>(dims_src));
    uint16_t* dst_ptr = dst.get();
    uint16 * dst_same_ptr = dst_same.get();
    uint16_t* src_ptr = src.get();
    uint16_t* empty_ptr = unique_ptr<uint16_t[]>().get();
    uint16_t dst_default_val = 13, src_default_val = 79;
    EXPECT_TRUE(!empty_ptr);
    EXPECT_TRUE(dst_ptr && src_ptr);
    vector<int> dims_wrong_size({1, 1, 1, 1}), bad_offsets({-1}),
                bad_extents({3, 0, 1}), bad_strides({2, 2, 0});

    // nullptr
    EXPECT_THROW(mcp3d::CopyDataVolume<uint16_t>(empty_ptr, dims_dst, src_ptr, dims_src),
                 mcp3d::MCPAssertionError);
    // bad dimension size
    EXPECT_THROW(mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_wrong_size, src_ptr, dims_src,
                                                 mcp3d::MImageBlock{},
                                                 mcp3d::MImageBlock({}, {}, strides_fit)),
                 mcp3d::MCPAssertionError);
    // bad offsets from MImageBlock
    EXPECT_THROW(mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                                 mcp3d::MImageBlock{},
                                                 mcp3d::MImageBlock(bad_offsets, dims_src, strides_unit)),
                 mcp3d::MCPOutOfRangeError);
    // dst too small to hold src. omitting all arguments with default value,
    // which copies src fully to dst with unit strides
    EXPECT_THROW(mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src),
                 mcp3d::MCPAssertionError);

    long addr_end;
    mcp3d::SetConstant(dst_ptr, dims_dst, dst_default_val);
    EXPECT_EQ(dst_default_val, dst_ptr[0]);
    mcp3d::SetConstant(dst_same_ptr, dims_dst_same, dst_default_val);
    EXPECT_EQ(dst_default_val, dst_same_ptr[0]);
    mcp3d::SetConstant(src_ptr, dims_src, src_default_val);
    EXPECT_EQ(src_default_val, src_ptr[0]);

    // copy continuous block
    mcp3d::CopyDataVolume<uint16_t>(dst_same_ptr, dims_dst_same,
                                    src_ptr, dims_src,
                                    mcp3d::MImageBlock{}, mcp3d::MImageBlock{});
    addr_end = mcp3d::ReduceProdSeq<long>(dims_dst_same);
    for (long i = 0; i < addr_end; ++i)
        EXPECT_EQ(src_default_val, dst_same_ptr[i]);

    // copy stride 1 block
    mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                    mcp3d::MImageBlock{},
                                    mcp3d::MImageBlock({5, 50, 50}, {5, 50, 50}));
    addr_end = mcp3d::ReduceProdSeq<long>(dims_dst);
    for (long i = 0; i < addr_end; ++i)
        EXPECT_EQ(src_default_val, dst_ptr[i]);

    // copy src data with exact fit stride
    mcp3d::SetConstant(dst_ptr, dims_dst, dst_default_val);
    mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                    mcp3d::MImageBlock{},
                                    mcp3d::MImageBlock({}, {}, strides_fit));
    addr_end = mcp3d::ReduceProdSeq<long>(dims_dst);
    for (long i = 0; i < addr_end; ++i)
        EXPECT_EQ(src_default_val, dst_ptr[i]);
    // copy single voxel, passing a larger than extents strides vector
    mcp3d::SetConstant(dst_ptr, dims_dst, dst_default_val);
    do
    {
        mcp3d::SetRandom<uint16_t>(src_ptr, dims_src);
    } while (src_ptr[0] == dst_ptr[0]);
    mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                    mcp3d::MImageBlock{},
                                    mcp3d::MImageBlock({0, 0, 0}, {1, 1, 1}, {2, 2, 2}));
    EXPECT_EQ(src_ptr[0], dst_ptr[0]);
    // copy from src to dst
    vector<int> src_offsets({2, 10, 10}), src_extents({4, 21, 21}), src_strides({1, 2, 2}),
                dst_offsets({1, 5, 5}), dst_strides({1, 4, 4});
    mcp3d::SetConstant(dst_ptr, dims_dst, dst_default_val);
    mcp3d::SetRandom<uint16_t>(src_ptr, dims_src);
    mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                    mcp3d::MImageBlock{dst_offsets, {}, dst_strides},
                                    mcp3d::MImageBlock(src_offsets, src_extents, src_strides));
    for (int z_dst = dst_offsets[0], z_src = src_offsets[0]; z_src < src_offsets[0] + src_extents[0]; z_dst += dst_strides[0], z_src += src_strides[0])
        for (int y_dst = dst_offsets[1], y_src = src_offsets[1]; y_src < src_offsets[1] + src_extents[1]; y_dst += dst_strides[1], y_src += src_strides[1])
        {
            int x_dst = dst_offsets[2], x_src = src_offsets[2];
            long addr_src = mcp3d::LinearAddress(dims_src, z_src, y_src, x_src),
                 addr_dst = mcp3d::LinearAddress(dims_dst, z_dst, y_dst, x_dst);
            for (; x_src < src_offsets[2] + src_extents[2]; x_src += src_strides[2])
            {
                EXPECT_EQ(src_ptr[addr_src], dst_ptr[addr_dst]);
                addr_dst += dst_strides[2];
                addr_src += src_strides[2];
            }
        }
    // should fail range check
    src_extents = {4, 23, 23};
    EXPECT_THROW(mcp3d::CopyDataVolume<uint16_t>(dst_ptr, dims_dst, src_ptr, dims_src,
                                                 mcp3d::MImageBlock{dst_offsets, {}, dst_strides},
                                                 mcp3d::MImageBlock(src_offsets, src_extents, src_strides)),
                 mcp3d::MCPAssertionError);

}

TEST(ImageUtils, DataVolumeEqual)
{
    vector<int> dims({1024, 1024, 3});
    size_t n = mcp3d::ReduceProdSeq<size_t>(dims);
    unique_ptr<double[]> ptr1 = make_unique<double[]>(n);
    unique_ptr<double[]> ptr2 = make_unique<double[]>(n);
    // ptr1 and ptr2 are identical
    mcp3d::SetRandom<double>(ptr1.get(), dims);
    memcpy(ptr2.get(), ptr1.get(), sizeof(double) * n);
    EXPECT_TRUE(mcp3d::VolumeEqual(ptr1.get(), ptr2.get(), dims, sizeof(double)));
    // ptr1 and ptr2 are different
    mcp3d::SetRandom<double>(ptr2.get(), dims);
    EXPECT_FALSE(
            mcp3d::VolumeEqual(ptr1.get(), ptr2.get(), dims, sizeof(double)));
}

TEST(ImageUtils, QuickSort)
{
    long n = 1000000;
    unique_ptr<int[]> input(new (nothrow) int[n]);
    MCP3D_ASSERT(input);
    mcp3d::MCPTensor1DMap<int> input_map(input.get(), n);
    input_map.setRandom();
    vector<int> input_copy;
    for (long i = 0; i < n; ++i)
        input_copy.push_back(input.get()[i]);
    INIT_TIMER(0)
    INIT_TIMER(1)
    TIC(0)
    mcp3d::QuickSort<int>(input.get(), n);
    TOC(0)
    TIC(1)
    sort(input_copy.begin(), input_copy.end());
    TOC(1)
    for (long i = 0; i < n; ++i)
        EXPECT_EQ(input_copy[i], input.get()[i]);
    cout << "time to array sort = " << ELAPSE(0) << "seconds" << endl;
    cout << "time to vector sort = " << ELAPSE(1) << "seconds" << endl;
}


