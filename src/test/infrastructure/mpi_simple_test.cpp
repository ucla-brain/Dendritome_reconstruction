#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <fstream>
#include <memory>
#include <mpi.h>
#include "common/mcp3d_common.hpp"
#include "common/mcp3d_utility.hpp"

using namespace std;

void write_mpi_rank_sum_output(int rank, int size,
                               string outfile_path, bool correct)
{
    char hostname[255];
    gethostname(hostname, 255);
    string result = correct? "process rank sum correct" : "process rank sum incorrect";
    ofstream f;
    f.exceptions(ofstream::failbit | ofstream::badbit);
    f.open(outfile_path, ios::app);
    f << "process " << to_string(rank) << '/' << to_string(size)
      << ": @" << string(hostname) << ", " << result << '\n';
    f.close();
}

bool mpi_rank_sum(int rank, int size, string outfile_path)
{
    bool correct, turn_to_write = false;
    int reduce_sum, correct_sum = 0, tag = 10;
    MPI_Status status;
    for (int i = 0; i < size; ++i)
        correct_sum += i;
    MPI_Allreduce(&rank, &reduce_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    correct = reduce_sum == correct_sum;
    if (rank == 0)
        MPI_Reduce(MPI_IN_PLACE, &correct, 1,
                   MPI_CXX_BOOL, MPI_LAND, 0, MPI_COMM_WORLD);
    else
        MPI_Reduce(&correct, &correct, 1,
                   MPI_CXX_BOOL, MPI_LAND, 0, MPI_COMM_WORLD);
    if (rank == 0)
        turn_to_write = true;
    // straight line send: 0 -> 1 -> ... rank_max
    // rank 0 does not receive, rank_max does not send
    // only a boolean sent, MPI_Send should return right away
    if (rank == 0)
    {
        write_mpi_rank_sum_output(rank, size, outfile_path, correct);
        MPI_Send(&turn_to_write, 1, MPI_CXX_BOOL, rank + 1, tag, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(&turn_to_write, 1, MPI_CXX_BOOL, rank - 1, tag, MPI_COMM_WORLD, &status);
        if (turn_to_write)
            write_mpi_rank_sum_output(rank, size, outfile_path, correct);
        if (rank < size - 1)
            MPI_Send(&turn_to_write, 1, MPI_CXX_BOOL, rank + 1, tag, MPI_COMM_WORLD);
    }
    return correct;
}

bool mpi_all_to_all(int rank, int size, string outfile_path)
{
    unique_ptr<int[]> send_buf(new int[size]);
    unique_ptr<int[]> rec_buf(new int[size]);
    for (int i = 0; i < size; ++i)
        send_buf[i] = rank;
    MPI_Alltoall(send_buf.get(), 1, MPI_INT, rec_buf.get(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 0; i < size; ++i)
        if (rec_buf[i] != i)
            return false;
    return true;
}

/*
 * test for MPI_Send, MPI_Recv, MPI_Reduce, MPI_Allreduce, and MPI_Alltoall
 * calls. the MPI program is aborted if results are incorrect.
 */
int main(int argc, char* argv[])
{
    int rank, size;
    string outfile_path = mcp3d::JoinPath(
            {mcp3d::test_data_dir(), "infrastructure",
             "mpi_simple_test.out"});

    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0 && mcp3d::IsFile(outfile_path))
        remove(outfile_path.c_str());
    MPI_Barrier(MPI_COMM_WORLD);

    bool rank_sum_correct = mpi_rank_sum(rank, size, outfile_path);
    if (!rank_sum_correct)
    {
        cout << "incorrect result in rank sum" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    bool all_to_all_correct = mpi_all_to_all(rank, size, outfile_path);
    if (rank == 0)
    {
        ofstream f;
        f.exceptions(ofstream::failbit | ofstream::badbit);
        f.open(outfile_path, ios::app);
        if (all_to_all_correct)
            f << "mpi all to all communication completes successfully" <<endl;
        else
            f << "mpi all to all communication error" <<endl;
    }
    if (!all_to_all_correct)
    {
        cout << "incorrect result in MPI all to all" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (rank == 0)
    {
        ifstream ifs(outfile_path);
        int count = 0;
        string line;
        while (getline(ifs, line))
            count++;
        if (count != size + 1)
        {
            cout << "expecting " << size + 1 << " lines of output" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Finalize();
}
