//
// Created by muyezhu on 10/8/17.
//
# define EIGEN_USE_THREADS
# define EIGEN_USE_BLAS
# define EIGEN_USE_LAPACKE
# define EIGEN_NO_DEBUG
# define EIGEN_UNALIGNED_VECTORIZE 1
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <mpi.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_intensity.hpp"

using namespace std;

static int process_id, n_proc;

string out_file_path(int id = -1)
{
    char host[255];
    gethostname(host, 255);
    if (id < 0)
        return mcp3d::JoinPath({mcp3d::benchmark_module_dir(),
                                "cluster_runtime", "cv_vs_eigen",
                                "cv_vs_eigen_" + string(host) + ".txt"});
    else
        return mcp3d::JoinPath({mcp3d::benchmark_module_dir(),
                                "cluster_runtime", "cv_vs_eigen",
                                "cv_vs_eigen_hpc_" + string(host) + "_size" +
                                to_string(n_proc) + "_rank" + to_string(id) +
                                ".txt"});
}

void opencv_2d_conv(const cv::Mat& m, const cv::Mat& g)
{
    cv::Mat mc;
    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    cv::sepFilter2D(m, mc, CV_64F, g, g);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    cout << "opencv 2d convolution: d0 = " << m.rows << ", d1 = " << m.cols
         << ", time = " << time_spent.count() << " milliseconds" << endl;
}

void eigen_2d_conv_default_device(const Eigen::Tensor<double, 2, Eigen::RowMajor>& t,
                                  const Eigen::Tensor<double, 1, Eigen::RowMajor>& kernel)
{
    Eigen::Tensor<double, 2, Eigen::RowMajor> tp, tc;
    Eigen::array<std::pair<int, int>, 2> paddings;
    int k = kernel.dimensions()[0];
    paddings[0] = std::make_pair(k / 2, k / 2);
    paddings[1] = std::make_pair(k / 2, k / 2);
    Eigen::array<int, 1> dim0 = {0};
    Eigen::array<int, 1> dim1 = {1};

    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    tp = t.pad(paddings);
    tc = tp.convolve(kernel, dim0)
           .convolve(kernel, dim1);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    cout << "eigen 2d convolution single thread: d0 = " << t.dimensions()[0] << ", d1 = "
         << t.dimensions()[1] << ", time = " << time_spent.count() << " milliseconds" << endl;
}

void eigen_2d_conv_simple_thread_pool(const Eigen::Tensor<double, 2, Eigen::RowMajor>& t,
                               const Eigen::Tensor<double, 1, Eigen::RowMajor>& kernel,
                               int n_thread)
{
    int k = kernel.dimensions()[0];
    Eigen::Tensor<double, 2, Eigen::RowMajor> tp (t.dimensions()[0] + k - 1, t.dimensions()[1] + k - 1);
    Eigen::Tensor<double, 2, Eigen::RowMajor> tc_0 (t.dimensions()[0], t.dimensions()[1] + k - 1);
    Eigen::Tensor<double, 2, Eigen::RowMajor> tc_1 (t.dimensions()[0], t.dimensions()[1]);
    Eigen::array<std::pair<int, int>, 2> paddings;
    paddings[0] = std::make_pair(k / 2, k / 2);
    paddings[1] = std::make_pair(k / 2, k / 2);
    Eigen::array<int, 1> dim0 = {0};
    Eigen::array<int, 1> dim1 = {1};

    Eigen::SimpleThreadPool interface(n_thread);
    Eigen::ThreadPoolDevice my_device(&interface, n_thread);
    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    tp.device(my_device) = t.pad(paddings);
    tc_0.device(my_device) = tp.convolve(kernel, dim0);
    tp.resize(0, 0);
    tc_1.device(my_device) = tc_0.convolve(kernel, dim1);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    cout << "eigen 2d convolution simple thread pool " << n_thread << " thread: d0 = " << t.dimensions()[0] << ", d1 = "
         << t.dimensions()[1] << ", time = " << time_spent.count() << " milliseconds" << endl;
}

void eigen_2d_conv_nonblock_thread_pool(const Eigen::Tensor<double, 2, Eigen::RowMajor>& t,
                                      const Eigen::Tensor<double, 1, Eigen::RowMajor>& kernel,
                                      int n_thread)
{
    int k = kernel.dimensions()[0];
    Eigen::Tensor<double, 2, Eigen::RowMajor> tp (t.dimensions()[0] + k - 1, t.dimensions()[1] + k - 1);
    Eigen::Tensor<double, 2, Eigen::RowMajor> tc_0 (t.dimensions()[0], t.dimensions()[1] + k - 1);
    Eigen::Tensor<double, 2, Eigen::RowMajor> tc_1 (t.dimensions()[0], t.dimensions()[1]);
    Eigen::array<std::pair<int, int>, 2> paddings;
    paddings[0] = std::make_pair(k / 2, k / 2);
    paddings[1] = std::make_pair(k / 2, k / 2);
    Eigen::array<int, 1> dim0 = {0};
    Eigen::array<int, 1> dim1 = {1};

    Eigen::NonBlockingThreadPool interface(n_thread);
    Eigen::ThreadPoolDevice my_device(&interface, n_thread);
    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    tp.device(my_device) = t.pad(paddings);
    tc_0.device(my_device) = tp.convolve(kernel, dim0);
    tp.resize(0, 0);
    tc_1.device(my_device) = tc_0.convolve(kernel, dim1);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    cout << "eigen 2d convolution non blocking thread pool " << n_thread << " thread: d0 = " << t.dimensions()[0] << ", d1 = "
         << t.dimensions()[1] << ", time = " << time_spent.count() << " milliseconds" << endl;
}

void opencv_3d_conv(const cv::Mat& m, int k, double sigma, ofstream &f)
{
    cv::Mat mc;
    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    cv::Mat gy = cv::getGaussianKernel(k, sigma);
    cv::Mat gx;
    cv::transpose(gy, gx);
    cv::sepFilter2D(m, mc, CV_32F, gy, gx);
    mc = mc.reshape(1, m.rows * m.cols);
    cv::filter2D(mc, mc, CV_32F, gx);
    mc = mc.reshape(m.channels(), m.rows);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    f << "opencv 3d convolution: " << time_spent.count() << " milliseconds" << endl;
}

// use 4 threads
void eigen_3d_conv(const mcp3d::MCPTensor3D<float> &data,
                   mcp3d::MCPTensor3D<float> &result,
                   int k, double sigma, ofstream &f, int n = 4)
{
    chrono::high_resolution_clock ::time_point tbegin = chrono::high_resolution_clock::now();
    mcp3d::gaussian_convolution_3d<float, float>(data, result, k, sigma, n);
    chrono::high_resolution_clock::time_point tend = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_spent = tend - tbegin;
    f << "eigen 3d convolution: " <<  time_spent.count() << " milliseconds" << endl;
}

void benchmark_conv2d_speed(int d0, int d1, int k, double sigma)

{
    Eigen::Tensor<double, 2, Eigen::RowMajor> t(d0, d1);
    t.setRandom<Eigen::internal::NormalRandomGenerator<double>>();
    cv::Mat m(d0, d1, CV_64F, t.data());

    cv::Mat g = cv::getGaussianKernel(k, sigma);
    Eigen::Tensor<double, 1, Eigen::RowMajor> kernel(k);
    for (int i = 0; i < k; ++i)
        kernel(i) = g.at<double>(i);
    cout << "array shape: (" << d0 << ", " << d1 << ")" << endl;
    cout << "k = " << k << endl;
    opencv_2d_conv(m, g);
    eigen_2d_conv_default_device(t, kernel);
    eigen_2d_conv_simple_thread_pool(t, kernel, 4);
    eigen_2d_conv_simple_thread_pool(t, kernel, 16);
    eigen_2d_conv_nonblock_thread_pool(t, kernel, 2);
    eigen_2d_conv_nonblock_thread_pool(t, kernel, 4);
    eigen_2d_conv_nonblock_thread_pool(t, kernel, 8);
    eigen_2d_conv_nonblock_thread_pool(t, kernel, 16);
    eigen_2d_conv_nonblock_thread_pool(t, kernel, 32);
}

void benchmark_conv3d_speed(int d0, int d1, int d2, int k, double sigma,
                            ofstream &f, int n = 4)

{
    mcp3d::MCPTensor3D<float> t(d0, d1, d2);
    t.setRandom<Eigen::internal::NormalRandomGenerator<float>>();
    mcp3d::MCPTensor3D<float> r;
    cv::Mat m(d0, d1, CV_32FC(d2), t.data());
    f << "array shape: (" << d0 << ", " << d1 << ", " << d2 << ")" << endl;
    f << "k = " << k << endl;
    f << "eigen device: non blocking thread pool with " << n << " threads" <<endl;
    opencv_3d_conv(m, k, sigma, f);
    eigen_3d_conv(t, r, k, sigma, f, n);
}

int main(int argc, char** argv)
{
    vector<int> ks({5, 15, 25, 27, 29, 31, 35, 45, 55});
    double sigma = 1.0;
    int n;
    if (argc == 3)
    {
        n = stoi(argv[1]);
        if (strcmp(argv[2], "parallel") == 0)
        {
            MPI_Init(nullptr, nullptr);
            MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
            MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
        }
        else
        {
            n_proc = 1;
            process_id = 0;
        }
    }

    else if (argc == 2)
    {
        n = stoi(argv[1]);
        n_proc = 1;
        process_id = 0;
    }

    else
    {
        n = 1;
        n_proc = 1;
        process_id = 0;
    }

    ofstream f;
    int f_id = (n_proc == 1) ? -1 : process_id;
    f.open(out_file_path(f_id), ios::out);
    MCP3D_ASSERT(f.good());
    char host[255];
    gethostname(host, 255);
    f << "@" << string(host);
    if (n_proc == 1)
        f << ": single core\n";
    else
        f << ": mpi process = " << n_proc << "\n";

    cout << cv::getBuildInformation() << endl;
    cout << "opencv parallel section thread number: " << cv::getNumThreads() << endl;
    for (const int k: ks)
    {
        benchmark_conv3d_speed(512, 512, k, k, sigma, f, n);
        benchmark_conv3d_speed(1024, 1024, k, k, sigma, f, n);
        benchmark_conv3d_speed(2048, 2048, k, k, sigma, f, n);
    }
    f.close();

    // clean up MPI environment
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized)
        MPI_Finalize();
    return 0;

}