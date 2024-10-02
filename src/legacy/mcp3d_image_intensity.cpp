//
// Created by muyezhu on 3/18/17.
//

#include <fstream>
#include <opencv2/imgproc/imgproc.hpp>
#include "common/mcp3d_utility.hpp"
#include "mcp3d_image_intensity.hpp"

using namespace std;

void mcp3d::gaussian_convolution_3d_dep(const cv::Mat &input,
                                        cv::Mat &output, int k, double sigma)
{
    int rows = input.rows;
    int cols = input.cols;
    int levels = input.channels();

    cv::Mat g = cv::getGaussianKernel(k, sigma, CV_64F);
    cv::Mat g_x;
    cv::transpose(g, g_x);

    cv::GaussianBlur(input, output, cv::Size(k, k), sigma);
    output.reshape(1, rows * cols);
    cv::filter2D(output, output, -1, g_x);
    output.reshape(levels, rows);
}

void mcp3d::non_linear_intensity_adj_3d_process(cv::Mat& stack, cv::Mat& C,
                                                double c1, double c2,
                                                int k, double sigma)
{
    double e = 1e-6;
    cv::Mat stack_convolv, B, T;

    double a1, a2;
    cv::minMaxLoc(stack, &a1, &a2);
    assert(a2 < 1 + e);

    // 3d convolution
    mcp3d::gaussian_convolution_3d_dep(stack, stack_convolv, k, sigma);

    // addition of epsilon to avoid division by zero
    cv::add(stack_convolv, e, stack_convolv);
    // brightness B = c1 * stack_convolv / (1 - stack_convolv) + c2
    cv::multiply(stack_convolv, -1, B);
    cv::add(B, 1, B);  // B = 1 - stack_convolv
    cv::divide(stack_convolv, B, B);  // B = stack_convolv / (1 - stack_convolv)
    cv::multiply(B, c1, B);
    cv::add(B, c2, B);   // B = c1 * stack_convolv / (1 - stack_convolv) + c2

    // contrast C = (stack / stack_convolv) * (T + T' * (stack - stack_convolv))
    // where T = power(stack_ave, B), T' = B * power(stack_ave, B - 1)
    mcp3d::mat_power(stack_convolv, B, T);  // T = power(stack_convolv, B)
    cv::Mat temp;
    cv::add(B, -1, temp);  // temp = B - 1
    mcp3d::mat_power(stack_convolv, temp, temp);  // temp = power(stack_convolv, B - 1)
    cv::multiply(B, temp, B);  // B = B * np.power(stack_convolv, B - 1)
    cv::subtract(stack, stack_convolv, temp);
    cv::multiply(B, temp, B); // B = B * np.power(stack_convolv, B - 1) * (stack - stack_convolv)
    cv::add(T, B, B);   // B = T + B * np.power(stack_convolv, B - 1) * (stack - stack_convolv)
    temp.release();
    T.release();
    cv::divide(stack, stack_convolv, C);  //C = stack / stack_convolv
    cv::multiply(C, B, C);   // C = (stack / stack_convolv) * (T + B * np.power(stack_convolv, B - 1) * (stack - stack_convolv))
    B.release();
    cv::min(C, 1, C);
    cv::max(C, 0, C);
}
/*
void mcp3d::non_linear_intensity_adj_3d(const vector<string>& img_paths,
                                        double c1, double c2,
                                        int k, double sigma)
{
    // load images in float64 datatype and normalized
    stack_loader::BSLCenteredAt loader(img_paths, k, CV_64F, true);
    vector<string> valid_paths = loader.get_valid_img_paths();

    uint32_t last_processed = 0;
    for (uint32_t i = 0; i < valid_paths.size(); i += k)
    {
        // load stack centered at valid_paths[i], with buffer images
        loader.build_stack_vec_centered_at(i);
        cv::Mat C;
        non_linear_intensity_adj_3d_process(loader.stack, C, c1, c2, k, sigma);

        pair<int, int> s_indices = loader.kernel_img_stack_levels();
        pair<int, int> l_indices = loader.kernel_img_list_indices();

        mcp3d::SaveStacksToImgs(C, valid_paths,
                                s_indices.first, s_indices.second,
                                l_indices.first, "adj", "intensity_adj/");
        last_processed = l_indices.second;
    }

    if (last_processed < img_paths.size() - 1)
    {
        loader.build_stack_vec_centered_at((int)img_paths.size() - 1);
        cv::Mat C;
        non_linear_intensity_adj_3d_process(loader.stack, C, c1, c2, k, sigma);

        pair<int, int> s_indices = loader.kernel_img_stack_levels();

        int last = (int)img_paths.size() - 1;
        mcp3d::SaveStacksToImgs(C, valid_paths,
                                s_indices.second - (last - last_processed - 1),
                                s_indices.second, last_processed + 1,
                                "adj", "intensity_adj/");
    }
}*/

void mcp3d::gaussian_convolution_3d(const mcp3d::Image3D& volume,
                                    mcp3d::Image3D& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    MCP3D_ASSERT(mcp3d::KnownVoxelType(volume.datatype()) &&
                 volume.data_ptr() != nullptr)
    // if datatype not initialized, initialize to uint8_t
    if (volume_convolve.datatype() == VoxelType::UNKNOWN)
        volume_convolve.Init(VoxelType::M8U, volume.dimensions());
    else MCP3D_ASSERT(volume.SameDimensions(volume_convolve))

    if (!mcp3d::UnsignedVoxeltype(volume.datatype()) &&
        !mcp3d::FloatingVoxeltype(volume.datatype())) MCP3D_DOMAIN_ERROR(
            "only supporting unsigned integer and floating point types")

    if (volume.datatype() == VoxelType::M8U)
        gaussian_convolution_3d<uint8_t, float>(volume.image<uint8_t>(),
                                              volume_convolve,
                                              k, sigma, n_threads);
    else if (volume.datatype() == VoxelType::M16U)
        gaussian_convolution_3d<uint16_t, float>(volume.image<uint16_t>(),
                                              volume_convolve,
                                              k, sigma, n_threads);
    else if (volume.datatype() == VoxelType::M32U)
        gaussian_convolution_3d<uint32_t, float>(volume.image<uint32_t>(),
                                              volume_convolve,
                                              k, sigma, n_threads);
    else if (volume.datatype() == VoxelType::M32F)
        gaussian_convolution_3d<float, float>(volume.image<float>(),
                                              volume_convolve,
                                              k, sigma, n_threads);
    else if (volume.datatype() == VoxelType::M64F)
        gaussian_convolution_3d<double, float>(volume.image<double>(),
                                               volume_convolve,
                                               k, sigma, n_threads);

    MCP3D_ASSERT(volume.SameDimensions(volume_convolve));
}



