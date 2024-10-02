//
// Created by mzhu on 10/30/17.
//
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <chrono>
#include <tuple>
#include <mpi.h>
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_intensity.hpp"

using namespace std;

static int process_id, n_proc;
const vector<int> ks({5, 15, 25, 35, 45, 55});
const vector<int> nthreads({1, 2, 4, 8});
const double sigma = 1.0;

typedef tuple<int, int, int> dim3;

string out_file_path(int id = -1)
{
    char host[255];
    gethostname(host, 255);
    if (id < 0)
        return mcp3d::JoinPath(
                {mcp3d::benchmark_module_dir(), "cluster_runtime",
                 "conv_nthreads_" + string(host) + ".txt"});
    else
        return mcp3d::JoinPath(
                {mcp3d::benchmark_module_dir(), "cluster_runtime",
                 "conv_nthreads_hpc" + string(host) + "_rank" + to_string(id) +
                 ".txt"});
}

void CombineOutput()
{
    string command = "cat " + mcp3d::benchmark_module_dir() +
                     "cluster_runtime/conv_nthreads_hpc*.txt > " +
                     mcp3d::benchmark_module_dir() +
                     "cluster_runtime/conv_nthreads_hpc.txt";
    system(command.c_str());
}

void FillDim3Vec(vector<dim3> &d3s)
{
    d3s.push_back(make_tuple(512, 512, 50));
    d3s.push_back(make_tuple(512, 512, 100));
    d3s.push_back(make_tuple(512, 512, 200));
    d3s.push_back(make_tuple(1024, 1024, 50));
    d3s.push_back(make_tuple(1024, 1024, 100));
    d3s.push_back(make_tuple(1024, 1024, 200));
    d3s.push_back(make_tuple(2048, 2048, 50));
    d3s.push_back(make_tuple(2048, 2048, 100));
    d3s.push_back(make_tuple(2048, 2048, 200));
}

void TimeConv(dim3 d3, int k, int n, ofstream &f)
{
    using namespace chrono;
    high_resolution_clock::time_point t_start, t_end;
    duration<double> t_delta;
    t_start = high_resolution_clock::now();
    mcp3d::MCPTensor3D<double> data(get<0>(d3), get<1>(d3), get<2>(d3));
    data.setRandom();
    t_end = high_resolution_clock::now();
    t_delta = duration_cast<duration<double>>(t_end - t_start);
    f << "initialize tensor: " << t_delta.count() << " seconds" << endl;
    mcp3d::MCPTensor3D<double> result;
    t_start = high_resolution_clock::now();
    mcp3d::gaussian_convolution_3d<double>(data, result, k, sigma, n);
    t_end = high_resolution_clock::now();
    t_delta = duration_cast<duration<double>>(t_end - t_start);
    f << "convolution with k = " << k << ", nthreads = " << n << ": " << t_delta.count() << " seconds" << endl;
}

int main(int argc, char** argv)
{
    // MPI_init if in parallel mode
    if (argc == 2 && strcmp(argv[1], "parallel") == 0)
    {
        MPI_Init(nullptr, nullptr);
        MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
        MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
        cout << "rank = " << process_id << ", size = " << n_proc << endl;
    }
    else
    {
        process_id = 0;
        n_proc = 1;
    }
    // set up tensor dimensions
    vector<dim3> d3s;
    FillDim3Vec(d3s);
    // set up output filestream, write hostname and process number
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
    // time convolution operations
    for (const dim3 &d3: d3s)
    {
        f << "tensor dimension: (" << get<0>(d3) << ", " << get<1>(d3) << ", " << get<2>(d3) << ')' << endl;
        for (const int k: ks)
            for (const int n: nthreads)
                TimeConv(d3, k, n, f);
    }
    f.close();
    // clean up MPI environment and combine outputs
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (process_id == 0)
            CombineOutput();
        MPI_Finalize();
    }

    return 0;
}
