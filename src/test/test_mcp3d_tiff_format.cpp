//
// Created by muyezhu on 2/26/18.
//
#include <memory>
#include <random>
#include <chrono>
#include <tiffio.h>
#include <boost/filesystem.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <gtest/gtest.h>
#include "common/mcp3d_common.hpp"
#include "image/mcp3d_channel_info.hpp"
#include "image/mcp3d_image.hpp"
#include "image/mcp3d_image_in_memory.hpp"
#include "image_interface/mcp3d_tiff_io.hpp"
#include "image/mcp3d_tiff_format.hpp"
#include "image/mcp3d_image_io.hpp"

using namespace std;

string test_tiff_format_dir(const string& subdir = string())
{
    string d = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                "image_formats", "tiff_format", subdir});
    MCP3D_ASSERT(mcp3d::IsDir(d))
    return d;
}

// under test_tiff_format_dir():
// Z0.tif, Z1.tif, Z2.tif. uint8 voxel type, 2048 x 2048 image dimensions
TEST(MTiffFormat, ReadInfoSingleLevel)
{
    mcp3d::MTiffFormat tiff_format;
    // read from directory
    vector<string> seq = mcp3d::FilesInDir(test_tiff_format_dir(), false, true,
                                           {".tif"});
    mcp3d::MChannelInfo channel_info = tiff_format.ReadChannelInfo(
            test_tiff_format_dir(), 0, true);

    EXPECT_EQ(size_t(1), channel_info.channel_pyr_infos().size());
    EXPECT_EQ(1, channel_info.n_pyr_infos());
    EXPECT_EQ(vector<int>({1}), channel_info.pyr_xy_ratios());
    EXPECT_EQ(test_tiff_format_dir(), channel_info.channel_root_dir());
    EXPECT_EQ(test_tiff_format_dir(), channel_info.channel_pyr_dir(0));
    EXPECT_EQ(seq, channel_info.channel_pyr_info(0).slice_image_names().at(""));
    EXPECT_EQ(mcp3d::FileFormat::TIFF,
              channel_info.channel_pyr_info(0).file_format());
    EXPECT_TRUE(channel_info.channel_pyr_info(0).is_level0());
    EXPECT_EQ(2048, channel_info.channel_pyr_info(0).xdim());
    EXPECT_EQ(2048, channel_info.channel_pyr_info(0).ydim());
    EXPECT_EQ(3, channel_info.channel_pyr_info(0).zdim());
    EXPECT_EQ(vector<int>({3, 2048, 2048}), channel_info.channel_pyr_info(0).xyz_dims());
    EXPECT_EQ(2048, channel_info.channel_pyr_info(0).chunk_xdim());
    EXPECT_EQ(2048, channel_info.channel_pyr_info(0).chunk_ydim());
    EXPECT_EQ(1, channel_info.channel_pyr_info(0).chunk_zdim());
    EXPECT_EQ(vector<int>({1, 2048, 2048}),
              channel_info.channel_pyr_info(0).chunk_xyz_dims());
    EXPECT_EQ(1, channel_info.channel_pyr_info(0).n_xchunks());
    EXPECT_EQ(1, channel_info.channel_pyr_info(0).n_ychunks());
    EXPECT_EQ(3, channel_info.channel_pyr_info(0).n_zchunks());
    EXPECT_EQ(3, channel_info.channel_pyr_info(0).n_total_chunks());
    EXPECT_EQ("zyx", channel_info.channel_pyr_info(0).dims_order());
    EXPECT_TRUE(channel_info.channel_pyr_info(0).is_level0());
    EXPECT_EQ(mcp3d::VoxelType::M8U, channel_info.channel_pyr_info(0).voxel_type());
}

// through MImageIO
// under JoinPath(test_tiff_format_dir(), chunk2048_uint8):
// 3 image levels with pyramid ratios 1, 2, 4
// level 0: Z0.tif, Z1.tif, Z2.tif. uint8 voxel type, 2048 x 2048 image dimensions
TEST(MTiffFormat, ReadInfoMultiLevels)
{
    string img_dir = test_tiff_format_dir("chunk2048_uint8");
    unique_ptr<mcp3d::MImageIO> io = make_unique<mcp3d::MImageIO>();
    vector<string> seq;
    mcp3d::MImageInfo image_info(img_dir, {""});
    mcp3d::MChannelInfo channel_info(io->ReadChannelInfo(image_info, 0));
    EXPECT_EQ(img_dir, channel_info.channel_root_dir());
    EXPECT_EQ(size_t(3), channel_info.channel_pyr_infos().size());
    EXPECT_EQ(3, channel_info.n_pyr_infos());
    EXPECT_EQ(vector<int>({1, 2 , 4}), channel_info.pyr_xy_ratios());
    for (int i = 0; i < channel_info.n_pyr_infos(); ++i)
    {
        int ratio = channel_info.pyr_xy_ratios()[i];
        EXPECT_EQ(mcp3d::FileFormat::TIFF, channel_info.file_format(i));
        EXPECT_EQ(2048 / ratio, channel_info.channel_pyr_infos()[i].xdim());
        EXPECT_EQ(2048 / ratio, channel_info.channel_pyr_infos()[i].ydim());
        EXPECT_EQ(3, channel_info.channel_pyr_infos()[i].zdim());
        EXPECT_EQ(vector<int>({3, 2048 / ratio, 2048 / ratio}), channel_info.xyz_dims(i));
        EXPECT_EQ(2048 / ratio, channel_info.channel_pyr_infos()[i].chunk_xdim());
        EXPECT_EQ(2048 / ratio, channel_info.channel_pyr_infos()[i].chunk_ydim());
        EXPECT_EQ(1, channel_info.channel_pyr_infos()[i].chunk_zdim());
        EXPECT_EQ(vector<int>({1, 2048 / ratio, 2048 / ratio}),
                  channel_info.channel_pyr_infos()[i].chunk_xyz_dims());
        EXPECT_EQ(1, channel_info.channel_pyr_infos()[i].n_xchunks());
        EXPECT_EQ(1, channel_info.channel_pyr_infos()[i].n_ychunks());
        EXPECT_EQ(3, channel_info.channel_pyr_infos()[i].n_zchunks());
        EXPECT_EQ(3, channel_info.channel_pyr_infos()[i].n_chunks());
        EXPECT_EQ("zyx", channel_info.channel_pyr_infos()[i].dims_order());
        EXPECT_EQ(i == 0, channel_info.channel_pyr_infos()[i].is_level0());
        EXPECT_EQ(mcp3d::VoxelType::M8U, channel_info.voxel_type(i));
        string pyr_level_dir;
        if (i == 0)
        {
            seq = mcp3d::FilesInDir(img_dir, false, true, {".tif"});
            pyr_level_dir = channel_info.channel_root_dir();
        }
        else
        {
            seq = mcp3d::FilesInDir(
                    mcp3d::JoinPath(img_dir, "pyr_level_" + mcp3d::PadNumStr(i, 2)),
                    false, true, {".tif"});
            pyr_level_dir = mcp3d::MChannelLayout::PyrLevelDirPath(
                    channel_info.channel_root_dir(), i);
        }
        EXPECT_EQ(pyr_level_dir, channel_info.channel_pyr_dir(i));
        EXPECT_EQ(seq, channel_info.channel_pyr_infos()[i].slice_image_names().at(""));
    }
}

// this only test for the strip / scan line based partial read write
void MTiffFormat_ReadPartialTiffImage_Impl(const vector<int> &dims,
                                           const vector<int> &offsets,
                                           const vector<int> &extents,
                                           const vector<int> &strides,
                                           bool is_tiled)
{
    // make a uint16 image of dims dimension, set its values to random, save to disk
    unique_ptr<uint16_t []> img (new uint16_t [mcp3d::ReduceProdSeq<size_t>(dims)]);
    uint16_t* img_ptr = img.get();
    string img_path = mcp3d::JoinPath(test_tiff_format_dir(), "partial_read",
                                      "random_" + to_string(dims[0]) + "_" +
                                      to_string(dims[1]) + ".tif");
    cv::Mat full_img;
    if (!is_tiled)
    {
        mcp3d::SetRandom<uint16_t>(img_ptr, dims);
        full_img = cv::Mat(dims[0], dims[1], CV_16U, img_ptr);
        cv::imwrite(img_path, full_img);
    }
    else
    {
        EXPECT_TRUE(dims[0] % 4 == 0);
        EXPECT_TRUE(dims[1] % 4 == 0);
        int tile_xdim = dims[1] / 4,
            tile_ydim = dims[0] / 4;
        unique_ptr<uint16_t []> tile_img =
                make_unique<uint16_t[]>((size_t)tile_xdim * (size_t)tile_ydim);

        ::TIFF* tif = ::TIFFOpen(img_path.c_str(), "w");
        mcp3d::SetTiledTiffTags(tif, (uint32_t)dims[1], (uint32_t)dims[0], (uint32_t)tile_ydim, (uint32_t)tile_xdim, (short)16, (short)1);
        for (int i = 0; i < 16; ++i)
        {
            mcp3d::SetRandom<uint16_t>(tile_img.get(), {1, tile_ydim, tile_xdim});
            ::TIFFWriteEncodedTile(tif, (uint32_t)i, tile_img.get(), tile_xdim * tile_ydim * sizeof(uint16_t));
        }
        ::TIFFClose(tif);
    }

    // read tiff info of the random image
    ::TIFF* tif = ::TIFFOpen(img_path.c_str(), "r");
    EXPECT_TRUE(tif);
    mcp3d::TiffDirectoryInfo tiff_info(tif);
    EXPECT_EQ(16, tiff_info.bits_per_sample);
    EXPECT_EQ(1, tiff_info.samples_per_pixel);
    ::TIFFClose(tif);

    // init timer
    INIT_TIMER(0)
    // read the whole tiff image
    TIC(0)
    full_img = cv::imread(img_path, cv::IMREAD_ANYDEPTH | cv::IMREAD_GRAYSCALE);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("time to cv::imread full image with dimensions " + mcp3d::JoinVector(dims, ", "), 0)

    // buffer for whole image
    uint16_t* full_img_ptr = full_img.ptr<uint16_t>() ;
    long addr_img, addr_full_image;
    int bv = 0;

    // MImage instance partial tiff read
    mcp3d::MImage mimage(mcp3d::ParentDir(img_path));
    mimage.ReadImageInfo(0);
    mimage.SelectView(mcp3d::MImageBlock(offsets, extents, strides), 0);
    EXPECT_TRUE(!mimage.image_info().empty());
    EXPECT_TRUE(dims == vector<int>({mimage.ydim(0), mimage.xdim(0)}));
    EXPECT_TRUE(mcp3d::StridedExtents(extents, strides) ==
                vector<int>({mimage.selected_view().ydim(0),
                             mimage.selected_view().xdim(0)}));

    // partial IO through mimage
    TIC(0)
    mimage.ReadData();
    TOC(0)
    uint16_t* mimage_ptr = mimage.Plane<uint16_t>(0, 0);
    EXPECT_TRUE(mimage_ptr);
    // check values
    for (int y_img = 0, y_full_img = offsets[0];
         y_img < mimage.selected_view().ydim(0);
         y_img += 1, y_full_img += strides[0])
    {
        int x_img = 0, x_full_img = offsets[1];
        addr_img = mcp3d::LinearAddress(mcp3d::StridedExtents(extents, strides),
                                        y_img, x_img);
        addr_full_image = mcp3d::LinearAddress(dims, y_full_img, x_full_img);
        for (; x_img < mimage.selected_view().xdim(0);
             x_img += 1, x_full_img += strides[1])
        {
            if (y_full_img >= dims[0] || x_full_img >= dims[1])
                EXPECT_EQ(bv, mimage_ptr[addr_img]);
            else
                EXPECT_EQ(full_img_ptr[addr_full_image], mimage_ptr[addr_img]);
            addr_img += 1;
            addr_full_image += strides[1];
        }
    }
    // remove random image
    EXPECT_TRUE(mcp3d::RemovePath(img_path));
    EXPECT_TRUE(!mcp3d::IsFile(img_path));
    REPORT_TIME_TO_COMPLETION("time to partial read through MImage: full image dimensions = " +
                              mcp3d::JoinVector(dims, ", ") + ", offsets = " +
                              mcp3d::JoinVector(offsets, ", ") + ", extents = " +
                              mcp3d::JoinVector(extents, ", ") + ", strides = " +
                              mcp3d::JoinVector(strides, ", "), 0)
}

TEST(MTiffFormat, ReadPartialTiffImage)
{
    vector<int> dims ({8192, 8192});
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> distribution_offset(0, 8191);
    uniform_int_distribution<int> distribution_extent(100, 2048);
    uniform_int_distribution<int> distribution_stride(1, 5);
    cout << "executing TiffFormat_ReadPartialTiffImage_Impl 5 times for stripped tiff image" << endl;
    for (int i = 0; i < 5; ++i)
    {
        // offsets and extents of the partial image view
        vector<int> offsets({distribution_offset(generator), distribution_offset(generator)}),
                    extents({distribution_extent(generator), distribution_extent(generator)}),
                    strides({distribution_stride(generator), distribution_stride(generator)});
        MTiffFormat_ReadPartialTiffImage_Impl(dims, offsets, extents, strides, false);
    }
    cout << "executing TiffFormat_ReadPartialTiffImage_Impl 5 times for tiled tiff image" << endl;
    for (int i = 0; i < 5; ++i)
    {
        // offsets and extents of the partial image view
        vector<int> offsets({distribution_offset(generator), distribution_offset(generator)}),
                    extents({distribution_extent(generator), distribution_extent(generator)}),
                    strides({distribution_stride(generator), distribution_stride(generator)});
        MTiffFormat_ReadPartialTiffImage_Impl(dims, offsets, extents, strides, true);
    }
}

// large_files directory has large rgb tiffs Z50.tif - Z55.tif
// read block: offsets {10000, 10000, 2}, extents {2048, 2048, 2}, strides {1, 1, 1}
// ground truth rgb image region saved as
// x10000_x2048_y10000_y2048/Z52_x10000_x2048_y10000_y2048.tif
// x10000_x2048_y10000_y2048/Z53_x10000_x2048_y10000_y2048.tif
// x10000_x2048_y10000_y2048/Z54_x10000_x2048_y10000_y2048.tif
TEST(MTiffFormat, ReadPartialTiffRGBImage)
{
    string img_root_dir(test_tiff_format_dir("large_files_read")),
           expected_dir("x10000_x2048_y10000_y2048");
    MCP3D_ASSERT(mcp3d::IsDir(mcp3d::JoinPath(img_root_dir, expected_dir)))
    mcp3d::MImage image(img_root_dir);
    image.ReadImageInfo(0);
    EXPECT_EQ(vector<int>({1, 1, 6, 26604, 22343}), image.dims());
    vector<int> offsets({2, 10000, 10000}),
                extents({3, 2048, 2048}),
                strides({1, 1, 1});
    mcp3d::MImageBlock image_block(offsets, extents, strides);
    image.SelectView(image_block, 0);
    image.ReadData();
    EXPECT_FALSE(image.loaded_view().empty());
    EXPECT_EQ(mcp3d::VoxelType::M8U, image.voxel_type());
    vector<string> loaded_basenames({"Z52", "Z53", "Z54"});
    unique_ptr<uint8_t[]> cv_ptr = make_unique<uint8_t []>(mcp3d::ReduceProdSeq<size_t>(extents));
    for (int i = 0; i < (int)loaded_basenames.size(); ++i)
    {
        int n_plane = extents[1] * extents[2];
        cv::Mat m = cv::imread(mcp3d::JoinPath(img_root_dir, expected_dir,
                                               loaded_basenames[i].append("_x10000_x2048_y10000_y2048.tif")),
                               cv::IMREAD_ANYDEPTH | cv::IMREAD_ANYCOLOR);
        EXPECT_EQ(CV_8UC3, m.type());
        // Note here CV_BGR2GRAY is used because opencv pack color bits in bgrbgr order,
        // different from rgbrgb pack order of libtiff
        #if CV_MAJOR_VERSION < 4
            cv::cvtColor(m, m, CV_BGR2GRAY);
        #else
            cv::cvtColor(m, m, cv::COLOR_BGR2GRAY);
        #endif
        EXPECT_EQ(CV_8U, m.type());
        memcpy(cv_ptr.get() + i * n_plane, m.ptr(), (size_t)n_plane);
    }
    mcp3d::MImageInMemory cv_wrapper;
    cv_wrapper.WrapData(vector<uint8_t*>({cv_ptr.get()}), extents);
    EXPECT_TRUE(image.HasEqualData(cv_wrapper));
}

TEST(MTiffFormat, WriteImagePyramid)
{
    string stripe_img_dir(test_tiff_format_dir("tile_and_stripe/stripe")),
           tile_img_dir(test_tiff_format_dir("tile_and_stripe/tile"));
    mcp3d::MImage img_stripe(stripe_img_dir), img_tile(tile_img_dir);
    img_stripe.ReadImageInfo(0);
    EXPECT_EQ(1, img_stripe.n_pyr_levels());
    img_stripe.WriteImagePyramid(0, 0, false, false, mcp3d::FileFormat::TIFF);
    EXPECT_EQ(2, img_stripe.n_pyr_levels());
    img_tile.ReadImageInfo(0);
    EXPECT_EQ(1, img_tile.n_pyr_levels());
    img_tile.WriteImagePyramid(0, 0, false, false, mcp3d::FileFormat::TIFF);
    EXPECT_EQ(2, img_tile.n_pyr_levels());
    img_stripe.SelectView(mcp3d::MImageBlock{}, 0, 1);
    img_stripe.ReadData();
    img_tile.SelectView(mcp3d::MImageBlock{}, 0, 1);
    img_tile.ReadData();
    EXPECT_TRUE(img_stripe.HasEqualData(img_tile));
    EXPECT_TRUE(mcp3d::RemovePath(mcp3d::JoinPath({stripe_img_dir, "pyr_level_01"})));
    EXPECT_TRUE(mcp3d::RemovePath(mcp3d::JoinPath({tile_img_dir, "pyr_level_01"})));

    // the img root dir has an image volume with compete pyramids
    string img_root_dir(test_tiff_format_dir("large_files_write"));
    // remove pyramid levels above 0
    vector<string> non_zero_pyr_level_dirs =
            mcp3d::MChannelLayout::NonZeroPyrLevelDirs(img_root_dir);
    for (const auto& pyr_level_dir: non_zero_pyr_level_dirs)
        EXPECT_TRUE(mcp3d::RemovePath(pyr_level_dir));

    mcp3d::MImage image(img_root_dir);
    image.ReadImageInfo(0);
    EXPECT_EQ(vector<int>({1, 1, 6, 26604, 22343}), image.dims());

    INIT_TIMER(0)

    int start_parent_level = 0,
        end_parent_level = 999;
    TIC(0)
    image.WriteImagePyramid(0, start_parent_level, false, true, mcp3d::FileFormat::TIFF);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("single threaded level 1 pyramid generation: \n"
                              "level 0 image dimensions =  " +
                              mcp3d::JoinVector(image.channel_pyr_info(0, 0).xyz_dims(), ", "), 0)
    TIC(0)
    image.WriteImagePyramids(0, start_parent_level + 1, end_parent_level, false,
                             true, mcp3d::FileFormat::TIFF);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("single threaded level 2 and all later levels pyramids generation: \n"
                              "level 0 image dimensions =  " +
                              mcp3d::JoinVector(image.channel_pyr_info(0, 0).xyz_dims(), ", "), 0)
    EXPECT_EQ(image.MaxNumPyramidLevels(0), image.n_pyr_levels(0));
    string channel_info_json_path(mcp3d::MChannelLayout::channel_info_path(image.channel_root_dir(0), 0));
    mcp3d::MChannelInfo saved_info(channel_info_json_path);
    EXPECT_EQ(saved_info, image.channel_info(0));

    TIC(0)
    image.WriteImagePyramid(0, 2, false, false, mcp3d::FileFormat::TIFF);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("multi-threading level 3 pyramid generation: level 2 image dimensions =  " +
                              mcp3d::JoinVector(image.channel_pyr_info(0, 0).xyz_dims(), ", "), 0)
    EXPECT_EQ(image.MaxNumPyramidLevels(0), image.n_pyr_levels(0));

    non_zero_pyr_level_dirs =
            mcp3d::MChannelLayout::NonZeroPyrLevelDirs(img_root_dir);
    for (const auto& pyr_level_dir: non_zero_pyr_level_dirs)
        EXPECT_TRUE(mcp3d::RemovePath(pyr_level_dir));
    EXPECT_TRUE(mcp3d::RemovePath(channel_info_json_path));
    image.image_info().ClearChannelInfos();
    image.ReadImageInfo(0);
    EXPECT_EQ(1, image.n_pyr_levels(0));
}



