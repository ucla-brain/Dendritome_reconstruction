//
// Created by muyezhu on 12/10/18.
//
#include <mpi.h>
#include "common/mcp3d_utility.hpp"

using namespace std;
int main(int argc, char** argv)
{
    MPI_Init(nullptr, nullptr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // handle invalid arguments in argv
    if (rank == 0)
    {
        if (argc < 2)
        {
            cout << "usage: romove_path_mpi root_dir" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    INIT_TIMER(0);
    string root_dir(argv[1]);
    if (mcp3d::IsFile(root_dir))
        mcp3d::RemovePath(root_dir);
    else
    {
        TIC_WORLD(0)
        mcp3d::AsyncRemovePathMPI(root_dir);
        TOC_WORLD(0)
    }
    if (rank == 0)
        REPORT_TIME_TO_COMPLETION("removing " + root_dir, 0)

    MPI_Finalize();
}
