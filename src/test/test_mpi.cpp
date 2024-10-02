//
// Created by muyezhu on 12/3/18.
//
#include <memory>
#include <chrono>
#include <thread>
#include <random>
#include <utility>
#include <unordered_set>
#include <mpi.h>
#include <gtest/gtest.h>
#include "image/mcp3d_image_utils.hpp"
#include "image/mcp3d_image.hpp"

using namespace std;


TEST(MPI, TestInit)
{
    int mpi;
    MPI_Initialized(&mpi);
    EXPECT_EQ(1, mpi);
}

// create communicator from ranks divisible by 2
// use comm_create_group such that the comm creation is only collective over
// qualifying ranks. this function will only finish for all ranks if the
// the communicator creation is not collective over disqualified ranks
TEST(MPI, CommCreateGroup)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int even_ranks[size / 2];
    for (int i = 0, j = 0; i < size; ++i)
        if (i % 2 == 0)
            even_ranks[j++] = i;

    if (rank % 2 == 0)
    {
        MPI_Group world_group, even_group;
        MPI_Comm even_comm;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group_incl(world_group, size / 2, even_ranks, &even_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, even_group, 0, &even_comm);
        int group_rank, group_size;
        MPI_Comm_rank(even_comm, &group_rank);
        MPI_Comm_size(even_comm, &group_size);
        EXPECT_EQ(size / 2, group_size);
        EXPECT_TRUE(group_rank <= rank);
        int rank_ = rank, group_rank_ = group_rank, rank_lower;
        MPI_Status status;
        if (group_rank == 0)
            MPI_Send(&group_rank_, 1, MPI_INT, group_rank + 1, 0, even_comm);
        else
        {
            MPI_Recv(&rank_lower, 1, MPI_INT, group_rank - 1, 0, even_comm, &status);
            EXPECT_EQ(rank_lower, group_rank - 1);
            if (group_rank < group_size - 1)
                MPI_Send(&group_rank_, 1, MPI_INT, group_rank + 1, 0, even_comm);
        }
        if (rank < size - 1)
            MPI_Send(&rank_, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    if (rank % 2 == 1)
    {
        int rank_lower;
        MPI_Status status;
        MPI_Recv(&rank_lower, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        EXPECT_EQ(rank - 1, rank_lower);
    }
}

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

pair<MPI_Comm, MPI_Group> CreateQualifiedComm(int n_disqualify, int* disqualified_ranks_array)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm qualified_comm;
    MPI_Group world_group, qualified_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_excl(world_group, n_disqualify, disqualified_ranks_array, &qualified_group);
    MPI_Comm_create_group(MPI_COMM_WORLD, qualified_group, 0, &qualified_comm);
    return make_pair(qualified_comm, qualified_group);
}

TEST(MPI, DiscardSlowest)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    this_thread::sleep_for (chrono::milliseconds(100 * rank));

    int n_disqualify = 1;
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
        pair<MPI_Comm, MPI_Group> qualified = CreateQualifiedComm(n_disqualify, disqualified_ranks_array);
        int size_qualified;
        MPI_Comm_size(qualified.first, &size_qualified);
        EXPECT_EQ(size - n_disqualify, size_qualified);
    }
    else
    {
        int disqualified_ranks_array[n_disqualify];
        bool qualify = true;
        MPI_Request request;
        SpeedTrial();
        MPI_Isend(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Bcast(disqualified_ranks_array, n_disqualify, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 0; i < n_disqualify; ++i)
            if (rank == disqualified_ranks_array[i])
                qualify = false;
        if (qualify)
        {
            pair<MPI_Comm, MPI_Group> qualified = CreateQualifiedComm(n_disqualify, disqualified_ranks_array);
            int size_qualified;
            MPI_Comm_size(qualified.first, &size_qualified);
            EXPECT_EQ(size - n_disqualify, size_qualified);
            EXPECT_NE(size - 1, rank);
        }
        else
        {
            MPI_Cancel(&request);
            EXPECT_EQ(size - 1, rank);
        }
    }
}


