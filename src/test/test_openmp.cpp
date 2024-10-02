//
// Created by muyezhu on 3/25/18.
//
#include <thread>
#include <chrono>
#include <omp.h>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_utility.hpp"
#include <gtest/gtest.h>

using namespace std;

TEST(OpenMP, Simple)
{
    CHECK_PARALLEL_MODEL
    cout << "hello world from single thread" << endl;
    uint32_t n_threads = thread::hardware_concurrency();
    vector<bool> thread_can_execute(n_threads, false);
    thread_can_execute[0] = true;
    #pragma omp parallel num_threads(n_threads)
    {
        int thread_id = omp_get_thread_num();
        if (thread_id > 0)
            EXPECT_FALSE(thread_can_execute[thread_id]);
        __omp_barrier__();
        std::this_thread::sleep_for (std::chrono::milliseconds(10 * thread_id));
        if (thread_can_execute[thread_id])
        {
            cout << "hello world from thread " << thread_id << endl;
            if ((uint32_t)thread_id < n_threads - 1)
            {
                thread_can_execute[thread_id + 1] = true;
                EXPECT_TRUE(thread_can_execute[thread_id + 1]);
            }
        }
    }
}

void OpenMPRank0CallAdd(vector<int>& values, int thread_id)
{
    EXPECT_EQ(thread::hardware_concurrency(), values.size());
    values[thread_id] += (thread_id + 1);
}

TEST(OpenMP, Rank0Call)
{
    CHECK_PARALLEL_MODEL
    uint32_t n_threads = thread::hardware_concurrency();
    vector<int> values(n_threads, 0), expected(n_threads, 0);
    expected[0] = 1;
    #pragma omp parallel num_threads(n_threads)
    {
        EXPECT_EQ(n_threads, omp_get_num_threads());
        int thread_id = omp_get_thread_num();
        RANK0_CALL(OpenMPRank0CallAdd, values, thread_id)
    }
    EXPECT_EQ(expected, values);
}

void OpenMPRank0CallSyncAdd(vector<int>& values, int thread_id)
{
    EXPECT_EQ(thread::hardware_concurrency(), values.size());
    if (thread_id == 0)
        this_thread::sleep_for(chrono::seconds(3));
    values[thread_id] += (thread_id + 1);
}

TEST(OpenMP, Rank0CallSync)
{
    CHECK_PARALLEL_MODEL
    uint32_t n_threads = thread::hardware_concurrency();
    vector<int> values(n_threads, 0), expected(n_threads, 0);
    expected[0] = 1;
    #pragma omp parallel num_threads(n_threads)
    {
        EXPECT_EQ(n_threads, omp_get_num_threads());
        int thread_id = omp_get_thread_num();
        RANK0_CALL_SYNC(OpenMPRank0CallSyncAdd, values, thread_id)
        if (values[0] == 0)
            values[thread_id] -= thread_id;
    }
    EXPECT_EQ(expected, values);
}

void ThrowingFunction(int thread_id)
{
    MCP3D_RUNTIME_ERROR(string("thread ").append(to_string(thread_id)).append(": throwing error"))
}

void OpenMPMultiThreadedExceptionsImpl()
{
    CHECK_PARALLEL_MODEL
    uint32_t n_threads = thread::hardware_concurrency();
    mcp3d::MultiThreadExceptions me;
    #pragma omp parallel num_threads(n_threads)
    {
        int thread_id = omp_get_thread_num();
        me.RunAndCaptureException(ThrowingFunction, thread_id);
        if (me.HasCapturedException())
        {
            #pragma omp cancel parallel
        }
        if (me.HasCapturedException())
            EXPECT_NO_FATAL_FAILURE("parallel region should have been canceled");
    }
    EXPECT_TRUE(me.HasCapturedException());
    MCP3D_RETHROW(me.e_ptr());
}

TEST(OpenMP, MultiThreadedExceptions)
{
    // execute 10 times
    for (int i = 0; i < 10; ++i)
        try
        {
            OpenMPMultiThreadedExceptionsImpl();
        }
        catch (...)
        {
            mcp3d::PrintNestedException(::std::current_exception(), __FILE__, __LINE__, __func__);
        }
}

