//
// Created by muyezhu on 11/29/18.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <memory>
#include <algorithm>
#include <tiffio.h>
#include <mpi.h>
#include "common/mcp3d_common.hpp"
#include "image/mcp3d_image.hpp"
#include "image/mcp3d_image_io.hpp"


using namespace std;

#define N_DISQUALIFY_CLUSTER 10
#define N_DISQUALIFY_LOCAL 1

void SpeedTrial()
{
    int n = 1024 * 1024 * 10;
    unique_ptr<double[]> data = make_unique<double[]>((size_t)(n));
    mcp3d::SetRandom<double>(data.get(), vector<int>({512, 512, 10}));
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    double* data_ptr = data.get();
    for (int i = 0; i < n; ++i)
    {
        double r = distribution(generator);
        data_ptr[i] *= r;
    }
}

MPI_Comm CreateQualifiedComm(int n_disqualify, int* disqualified_ranks_array)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm qualified_comm;
    MPI_Group world_group, qualified_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_excl(world_group, n_disqualify, disqualified_ranks_array, &qualified_group);
    MPI_Comm_create_group(MPI_COMM_WORLD, qualified_group, 0, &qualified_comm);
    return qualified_comm;
}

int main(int argc, char** argv)
{
    MPI_Init(nullptr, nullptr);
    ::TIFFSetWarningHandler(nullptr);
    try
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        string home_dir(getenv("HOME"));
        time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        char time_cstr[13];
        strftime(time_cstr, 13, "%Y%m%d%H%M", localtime(&now));
        string err_log_path = mcp3d::JoinPath({home_dir, string("write_image_pyramid_err") + time_cstr});
        ofstream ofs;
        // handle invalid arguments in argv, logging to HOME directory
        if (rank == 0)
        {
            ofs.open(err_log_path, ofstream::out | ofstream::app);
            if (argc < 4)
            {
                ofs << "usage: write_image_pyramid img_root_dir parent_level output_format [abort_all_on_failure]" << endl;
                ofs.close();
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        string img_root_dir(argv[1]);
        if (rank == 0)
        {
            if (!mcp3d::IsDir(img_root_dir))
            {
                ofs.open(err_log_path, ofstream::out | ofstream::app);
                ofs << img_root_dir + " is not a directory" << endl;
                ofs.close();
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        int parent_level = 0;
        try
        {
            parent_level = stoi(string(argv[2]));
        }
        catch(const invalid_argument& e)
        {
            if (rank == 0)
            {
                ofs.open(err_log_path, ofstream::out | ofstream::app);
                ofs << "can not parse argument parent_level to integer" << endl;
                ofs.close();
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        string out_format_str(argv[3]);
        mcp3d::FileFormat out_format = mcp3d::FileFormatExtToEnum(
                out_format_str);
        bool abort_all_on_fail = argc >= 5 ?
                                 mcp3d::StringLower(string(argv[4])) == "true" : true;

        MPI_Barrier(MPI_COMM_WORLD);

        // try to discard very slow workers
        int n_disqualify = mcp3d::HostOnCluster() ? N_DISQUALIFY_CLUSTER : N_DISQUALIFY_LOCAL;
        if (rank == 0)
        {
            int sender_rank;
            MPI_Status status;
            unordered_set<int> disqualified_ranks;
            int disqualified_ranks_array[n_disqualify];
            for (int i = 1; i < size; ++i)
                disqualified_ranks.insert(i);
            while (disqualified_ranks.size() > (size_t)n_disqualify)
            {
                MPI_Recv(&sender_rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                disqualified_ranks.erase(sender_rank);
            }
            int i = 0;
            for (const auto& disqualify_rank: disqualified_ranks)
                disqualified_ranks_array[i++] = disqualify_rank;
            MPI_Bcast(disqualified_ranks_array, n_disqualify, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Comm qualified_comm = CreateQualifiedComm(n_disqualify, disqualified_ranks_array);
            unique_ptr<mcp3d::MImageIO> image_io = make_unique<mcp3d::MImageIO>();
            image_io->WriteImagePyramidMPI(img_root_dir, 0, parent_level,
                                           out_format, abort_all_on_fail,
                                           qualified_comm);
        }
        else
        {
            int disqualified_ranks_array[n_disqualify];
            bool qualify = true;
            MPI_Request request;
            SpeedTrial();
            MPI_Isend(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            MPI_Bcast(disqualified_ranks_array, n_disqualify, MPI_INT, 0,
                      MPI_COMM_WORLD);
            for (int i = 0; i < n_disqualify; ++i)
                if (rank == disqualified_ranks_array[i])
                    qualify = false;
            if (qualify)
            {
                MPI_Comm qualified_comm = CreateQualifiedComm(n_disqualify,
                                                              disqualified_ranks_array);
                unique_ptr<mcp3d::MImageIO> image_io = make_unique<mcp3d::MImageIO>();
                image_io->WriteImagePyramidMPI(img_root_dir, 0, parent_level,
                                               out_format, abort_all_on_fail,
                                               qualified_comm);
            }
            else
                MPI_Cancel(&request);
        }
    }
    catch (...)
    {
        MCP3D_PRINT_NESTED_EXCEPTION
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

