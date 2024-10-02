//
// Created by muyezhu on 3/18/17.
//

#ifndef MCP3D_IMG_PROCESSING_HPP
#define MCP3D_IMG_PROCESSING_HPP

#include <type_traits>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <unsupported/Eigen/CXX11/ThreadPool>
#include "common/mcp3d_common.hpp"
#include "mcp3d_image_common.hpp"
#include "mcp3d_image.hpp"
#include "pixel_conversion.hpp"

namespace mcp3d
{

/// a non linear adaptive method to adjust brightness and contrast.
void non_linear_intensity_adj_3d(const std::vector<std::string>& img_paths,
                                 double c1, double c2,
                                 int k = 21, double sigma = 1.0);

void non_linear_intensity_adj_3d_process(cv::Mat& stack, cv::Mat& C,
                                         double c1, double c2, int k = 21,
                                         double sigma = 1.0);

/// cv::Mat container based 3d convolution. depcreated. only used in funcitions
/// written before gaussian_convolution_3d
void gaussian_convolution_3d_dep(const cv::Mat &input, cv::Mat &output,
                                 int k, double sigma);

/// Image3D -> convolve to Image3D.
/// currently volume_convolve is float32 datatype
void gaussian_convolution_3d(const Image3D& volume,
                             Image3D& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

template<typename T1, typename T2>
void gaussian_convolution_3d(const MCPTensor3DMap<T1>& volume,
                             Image3D& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

template<typename T1, typename T2>
void gaussian_convolution_3d(const MCPTensor3D<T1>& volume,
                             Image3D& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

template<typename T1, typename T2>
void gaussian_convolution_3d(const MCPTensor3DMap<T1>& volume,
                             MCPTensor3D<T2>& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

/// if output to Tensor, Tensor datatype must be floating point
template<typename T1, typename T2>
void gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T1>& volume,
                             mcp3d::MCPTensor3DMap<T2>& volume_convolve,
                             int k, double sigma, int n_threads);

template<typename T>
void gaussian_convolution_3d(const MCPTensor3DMap<T>& volume,
                             MCPTensor3DMap<double>& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

template<typename T>
void gaussian_convolution_3d(const MCPTensor3DMap<T>& volume,
                             MCPTensor3DMap<float>& volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

template<typename T1, typename T2>
void gaussian_convolution_3d(const MCPTensor3D<T1> &volume,
                             MCPTensor3D<T2> &volume_convolve,
                             int k, double sigma,
                             int n_threads = N_THREADS);

}


template<typename T1, typename T2>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T1>& volume,
                                    mcp3d::MCPTensor3DMap<T2>& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    if (std::is_same<float, T2>:: value || std::is_same<double, T2>:: value)
        mcp3d::gaussian_convolution_3d(volume, volume_convolve, k, sigma, n_threads);
    else MCP3D_DOMAIN_ERROR("convolution output should be floating point")
}

template<typename T>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T>& volume,
                                    mcp3d::MCPTensor3DMap<float>& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    if (k % 2 == 0) MCP3D_INVALID_ARGUMENT("k must be an odd number")
    if (k < 3) MCP3D_INVALID_ARGUMENT("k must be at least 3")

    MCP3D_ASSERT(mcp3d::SameDimensions(volume, volume_convolve))

    // non blocking thread pool
    Eigen::NonBlockingThreadPool interface(n_threads);
    Eigen::ThreadPoolDevice thread_pool(&interface, n_threads);
    Eigen::array<std::pair<int, int>, 3> paddings;
    paddings[0] = std::make_pair(k / 2, k / 2);
    paddings[1] = std::make_pair(k / 2, k / 2);
    paddings[2] = std::make_pair(k / 2, k / 2);
    Eigen::array<int, 3> shrink_dimension = {1, 1, 1};
    mcp3d::MCPTensor3D<float> volume_pad(volume.dimensions()[0] + k - 1,
                                         volume.dimensions()[1] + k - 1,
                                         volume.dimensions()[2] + k - 1);
    if (! std::is_same<T, float>::value)
    {
        mcp3d::MCPTensor3D<float> volume_cast(volume.dimensions()[0],
                                              volume.dimensions()[1],
                                              volume.dimensions()[2]);
        volume_cast.device(thread_pool) = volume.template cast<float>();
        volume_pad.device(thread_pool) = volume_cast.pad(paddings);
        volume_cast.resize(shrink_dimension);
    }

    else
        volume_pad.device(thread_pool) = volume.template cast<float>()
                                               .pad(paddings);

    // set up kernel and dimensions
    Eigen::Tensor<float, 2, Eigen::RowMajor> kernel(k, 1);
    cv::Mat gaussian1d = cv::getGaussianKernel(k, sigma, CV_32F);
    for (int i = 0; i < k; ++i)
        kernel(i) = gaussian1d.at<float>(i);
    Eigen::array<int, 1> dim0 = {0};
    Eigen::array<int, 1> dim1 = {1};
    Eigen::array<int, 1> dim2 = {2};

    // convolve dim2
    mcp3d::MCPTensor3D<float> volume_c2(volume.dimensions()[0] + k - 1,
                                        volume.dimensions()[1] + k - 1,
                                        volume.dimensions()[2]);
    volume_c2.device(thread_pool) = volume_pad.convolve(kernel, dim2);
    volume_pad.resize(volume.dimensions()[0] + k - 1,
                      volume.dimensions()[1],
                      volume.dimensions()[2]);

    // convolve dim1
    volume_pad.device(thread_pool) = volume_c2.convolve(kernel, dim1);
    volume_c2.resize(shrink_dimension);

    // convolve dim0
    volume_convolve.device(thread_pool) = volume_pad.convolve(kernel, dim0);
    volume_pad.resize(shrink_dimension);

    MCP3D_ASSERT(mcp3d::SameDimensions(volume, volume_convolve))
}

template<typename T>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T>& volume,
                                 mcp3d::MCPTensor3DMap<double>& volume_convolve,
                                 int k, double sigma, int n_threads)
{
    if (k % 2 == 0) MCP3D_INVALID_ARGUMENT("k must be an odd number")
    if (k < 3) MCP3D_INVALID_ARGUMENT("k must be at least 3")

    MCP3D_ASSERT(mcp3d::SameDimensions(volume, volume_convolve))

    // non blocking thread pool
    Eigen::NonBlockingThreadPool interface(n_threads);
    Eigen::ThreadPoolDevice thread_pool(&interface, n_threads);
    Eigen::array<std::pair<int, int>, 3> paddings;
    paddings[0] = std::make_pair(k / 2, k / 2);
    paddings[1] = std::make_pair(k / 2, k / 2);
    paddings[2] = std::make_pair(k / 2, k / 2);
    Eigen::array<int, 3> shrink_dimension = {1, 1, 1};
    mcp3d::MCPTensor3D<double> volume_pad(volume.dimensions()[0] + k - 1,
                                          volume.dimensions()[1] + k - 1,
                                          volume.dimensions()[2] + k - 1);
    // pad input. volume is dependent name
    volume_pad.device(thread_pool) = volume.template cast<double>()
                                           .pad(paddings);
    // set up kernel and dimensions
    Eigen::Tensor<double, 2, Eigen::RowMajor> kernel(k, 1);
    cv::Mat gaussian1d = cv::getGaussianKernel(k, sigma, CV_64F);
    for (int i = 0; i < k; ++i)
        kernel(i) = gaussian1d.at<double>(i);
    Eigen::array<int, 1> dim0 = {0};
    Eigen::array<int, 1> dim1 = {1};
    Eigen::array<int, 1> dim2 = {2};

    // convolve dim2
    mcp3d::MCPTensor3D<double> volume_c2(volume.dimensions()[0] + k - 1,
                                         volume.dimensions()[1] + k - 1,
                                         volume.dimensions()[2]);
    volume_c2.device(thread_pool) = volume_pad.convolve(kernel, dim2);
    volume_pad.resize(volume.dimensions()[0] + k - 1,
                      volume.dimensions()[1],
                      volume.dimensions()[2]);

    // convolve dim1
    volume_pad.device(thread_pool) = volume_c2.convolve(kernel, dim1);
    volume_c2.resize(shrink_dimension);

    // convolve dim0
    volume_convolve.device(thread_pool) = volume_pad.convolve(kernel, dim0);
    volume_pad.resize(shrink_dimension);

    MCP3D_ASSERT(mcp3d::SameDimensions(volume, volume_convolve))
}

template<typename T1, typename T2>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T1>& volume,
                                    mcp3d::MCPTensor3D<T2>& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    if (k % 2 == 0) MCP3D_INVALID_ARGUMENT("k must be an odd number")
    if (k < 3) MCP3D_INVALID_ARGUMENT("k must be at least 3")
    // volume_convolve is either size 0 or same dimensions as volume
    if (volume_convolve.dimensions()[0] == 0)
        volume_convolve.resize(volume.dimensions()[0],
                               volume.dimensions()[1],
                               volume.dimensions()[2]);
    MCP3D_ASSERT(mcp3d::SameDimensions(volume, volume_convolve))
    mcp3d::MCPTensor3DMap<T2> volume_convolve_map(volume_convolve.data(),
                                                  volume.dimensions()[0],
                                                  volume.dimensions()[1],
                                                  volume.dimensions()[2]);
    mcp3d::gaussian_convolution_3d<T1, T2>(volume, volume_convolve_map,
                                           k, sigma, n_threads);
}

template<typename T1, typename T2>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3DMap<T1>& volume,
                                    mcp3d::Image3D& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    if (k % 2 == 0) MCP3D_INVALID_ARGUMENT("k must be an odd number")
    if (k < 3) MCP3D_INVALID_ARGUMENT("k must be at least 3")
    if (!mcp3d::KnownVoxelType(volume_convolve.datatype())) MCP3D_INVALID_ARGUMENT(
            "volume_convolve must have initialized data type")
    if (! volume_convolve.SameDimensions(volume)) MCP3D_INVALID_ARGUMENT(
            "volume_convolve must have same dimensions as volume")

    mcp3d::MCPTensor3D<T2> volume_convolve_float(volume.dimensions()[0],
                                                 volume.dimensions()[1],
                                                 volume.dimensions()[2]);
    mcp3d::gaussian_convolution_3d<T1, T2>(volume, volume_convolve_float,
                                           k, sigma, n_threads);
    mcp3d::saturation_cast(volume_convolve_float, volume_convolve);
}

template<typename T1, typename T2>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3D<T1>& volume,
                                    mcp3d::Image3D& volume_convolve,
                                    int k, double sigma, int n_threads)
{
    MCP3D_ASSERT(volume.data())
    mcp3d::MCPTensor3DMap<T1> volume_map((T1*)volume.data(),
                                         volume.dimensions()[0],
                                         volume.dimensions()[1],
                                         volume.dimensions()[2]);
    mcp3d::gaussian_convolution_3d<T1, T2>(volume_map, volume_convolve,
                                           k, sigma, n_threads);
}

template<typename T1, typename T2>
void mcp3d::gaussian_convolution_3d(const mcp3d::MCPTensor3D<T1> &volume,
                                    mcp3d::MCPTensor3D<T2> &volume_convolve,
                                    int k, double sigma,
                                    int n_threads)
{
    MCP3D_ASSERT(volume.data())
    mcp3d::MCPTensor3DMap<T1> volume_map((T1*)volume.data(),
                                        volume.dimensions()[0],
                                        volume.dimensions()[1],
                                        volume.dimensions()[2]);
    mcp3d::gaussian_convolution_3d<T1, T2>(volume_map, volume_convolve,
                                           k, sigma, n_threads);
}

#endif //MCP3D_IMG_PROCESSING_HPP
