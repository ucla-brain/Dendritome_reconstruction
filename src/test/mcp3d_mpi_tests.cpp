//
// Created by muyezhu on 12/3/18.
//
#include <iostream>
#include <mpi.h>
#include <gtest/gtest.h>

using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(nullptr, nullptr);
    ::testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    return MPI_Finalize();
}
