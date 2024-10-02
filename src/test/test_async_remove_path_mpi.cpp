//
// Created by muyezhu on 12/10/18.
//
#include <thread>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include "common/mcp3d_paths.hpp"
#include "common/mcp3d_utility.hpp"
#include "test_utility_internal.hpp"

using namespace std;

TEST(MPIOperation, AsyncRemovePath)
{
    // create random directory structures
    string async_remove_dir = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                             "utility", "async_remove");
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MCP3D_MESSAGE("rank = " + to_string(rank))

    int n_items;
    long seed;
    INIT_TIMER(0)
    cout << "creating 5 random directory structures for multiprocess deletion. n processor = " << size << endl;
    for (int i = 0; i < 5; ++i)
    {
        seed = chrono::system_clock::now().time_since_epoch().count();
        mcp3d::MakeDirectories(async_remove_dir);
        EXPECT_TRUE(boost::filesystem::is_directory(async_remove_dir));
        TIC_WORLD(0)
        if (rank == 0)
        {
            n_items = mcp3d::test::CreateDirectoryTree(async_remove_dir, seed);
            MPI_Bcast(&n_items, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        else
            MPI_Bcast(&n_items, 1, MPI_INT, 0, MPI_COMM_WORLD);
        TOC_WORLD(0)
        cout << "directory " << async_remove_dir << " contains " << n_items << " items" << endl;
        REPORT_TIME_TO_COMPLETION("time to create directory structure", 0)
        TIC_WORLD(0)
        mcp3d::AsyncRemovePathMPI(async_remove_dir);
        TOC_WORLD(0)
        REPORT_TIME_TO_COMPLETION("time to remove directory structure with " + to_string(size) + " processors", 0)
        EXPECT_FALSE(boost::filesystem::is_directory(async_remove_dir));
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

