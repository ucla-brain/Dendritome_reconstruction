//
// Created by muyezhu on 1/19/18.
//
#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <vector>
#include <unordered_map>
#include <mpi.h>
#include <boost/filesystem.hpp>
#include "common/mcp3d_common.hpp"
#include "parallel/mpio_util.hpp"
#include "common/mcp3d_utility.hpp"

using namespace std;
using namespace chrono;

void IndividualProcessWrite(char* buffer, streamsize n,
                            const string &writer_path, fstream &logger)
{
    fstream writer(writer_path.c_str(), fstream::out);
    auto start = high_resolution_clock::now();
    writer.write(buffer, n);
    writer.close();
    auto end = high_resolution_clock::now();
    duration<double, milli> elapse = end - start;
    logger << "individual process writing outside of MPI:"
           << elapse.count() << "ms" << endl;
}

int main(int argc, char **argv)
{
    if (argc != 2) MCP3D_INVALID_ARGUMENT(
            "usage: ./mpi_io_test romio_json_path")
    string romio_json_path(argv[1]);
    MPI_Init(nullptr, nullptr);
    int rank, size, writers_rank, writers_size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // make directories
    string out_dir = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                      "infrastructure", "mpi_io",
                                      to_string(size) + "writers",
                                      "output"});
    string out_log_dir = mcp3d::JoinPath({out_dir, "logs"});
    if (rank == 0)
    {
        using namespace boost::filesystem;
        if (mcp3d::IsDir(out_dir))
            remove_all(path(out_dir));
        mcp3d::MakeDirectories(out_dir);
        mcp3d::MakeDirectories(out_log_dir);
        if (!is_directory(path(out_dir)) || !is_directory(path(out_log_dir)))
        {
            cout << "directories not created. aborting." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // create communicator for writers
    MPI_Comm writers;
    int color = (rank % 1 == 0) ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &writers);
    if (writers != MPI_COMM_NULL)
    {
        MPI_Comm_size(writers, &writers_size);
        MPI_Comm_rank(writers, &writers_rank);
    }

    // set up output files
    string out_file_base = mcp3d::JoinPath({out_dir, "mpi_io_ints_"});
    string individual_file_path =
            mcp3d::JoinPath({out_dir, to_string(rank) + ".out"});
    string process_log_path = mcp3d::JoinPath({out_log_dir,
                                               "mpi_io_ints.process" +
                                               to_string(rank) + "_log"});
    string master_log_path = mcp3d::JoinPath({out_log_dir,
                                              "mpi_io_ints.master_log"});
    fstream process_log, master_log;
    process_log.open(process_log_path.c_str(), fstream::out);

    // MPI types
    MPI_Offset offset;
    MPI_Status status;
    MPI_Request request;
    // MPI info objects
    MPI_Info info;
    mcp3d::SetMPIInfo(info, romio_json_path);

    high_resolution_clock::time_point start, end;
    duration<double, milli> elapse;
    vector<long> length_vec({1000l, 1000000l, 100000000l});
    long max_per_write_length = 100000l;
    // master logging
    if (rank == 0)
    {
        master_log.open(master_log_path.c_str(), fstream::out);
        mcp3d::LogMPIInfo(info, master_log);
        master_log << "number of process = " << size
                   << ", number of writers = " << writers_size << endl;
        master_log << "\nnon collective writes" << endl;
    }
    // process logging
    mcp3d::LogMPIInfo(info, process_log);
    if (writers != MPI_COMM_NULL)
    {
        process_log << "number of process = " << size
                    << ", writing process = " << writers_size << endl;
        process_log << "process rank = " << rank << ", host = "
                    << mcp3d::HostName() << endl;
        process_log << "\nnon collective writes" << endl;
    }
    else
    {
        process_log << "process rank = " << rank << ", not a writer, host = "
                    << mcp3d::HostName() << endl;
        cout << "not a writer, rank = " << rank << endl;
        MPI_Finalize();
        exit(0);
    }

    for (const long &length: length_vec)
    {
        MPI_File file;
        unique_ptr<int32_t[]> buffer(new int32_t[length]);
        // each process write n = length number of integer equal to its rank
        for (int i = 0; i < length; ++i)
            buffer[i] = rank;
        process_log << "writing " << length << " integers" << endl;
        MPI_Barrier(writers);
        start = high_resolution_clock::now();
        IndividualProcessWrite((char*)(buffer.get()), length * sizeof(int32_t),
                               individual_file_path, process_log);
        MPI_Barrier(writers);
        end = high_resolution_clock::now();
        elapse = end - start;
        process_log << "individual processes writing outside of MPI:"
                    << elapse.count() << "ms" << endl;
        if (rank == 0)
            master_log << "writing individual file outside of MPI:"
                       << elapse.count() << "ms" << endl;
        boost::filesystem::remove(boost::filesystem::path(individual_file_path));
        string out_file_non_collective(out_file_base);
        out_file_non_collective.append(to_string(length)).append("_non_collective");
        if (mcp3d::IsFile(out_file_non_collective))
            remove(out_file_non_collective.c_str());
        // create file
        int error = MPI_File_open(
                writers, out_file_non_collective.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
                info, &file);
        if (error)
            MPI_Abort(writers, mcp3d::MPI_FILE_CREATE_ERR);
        // MPI array type for file view
        MPI_Datatype arraytype;
        MPI_Type_contiguous(length, MPI_INT32_T, &arraytype);
        MPI_Type_commit(&arraytype);
        offset = long(writers_rank) * length * sizeof(int32_t);
        MPI_File_set_view(file, offset, MPI_INT32_T, arraytype,
                          "native", info);
        MPI_File_set_atomicity(file, false);
        // MPI write
        MPI_Barrier(writers);
        start = high_resolution_clock::now();
        MPI_File_iwrite(file, buffer.get(), length, MPI_INT32_T, &request);
        MPI_Wait(&request, &status);
        end = high_resolution_clock::now();
        elapse = end - start;
        process_log << "file write:" << elapse.count() << "ms" << endl;
        MPI_Barrier(writers);
        end = high_resolution_clock::now();
        elapse = end - start;
        if (rank == 0)
            master_log << "writing " << length << " integers, time = "
                       << elapse.count() << " milliseconds\n" << endl;
        MPI_File_close(&file);
    }

    // collective
    if (rank == 0)
        master_log << "\ncollective writes" << endl;
    process_log << "\ncollective writes" << endl;
    for (const long &length: length_vec)
    {
        MPI_File file;
        unique_ptr<int32_t[]> buffer(new int32_t[length]);
        // each process write n = length number of integer equal to its rank
        for (int i = 0; i < length; ++i)
            buffer[i] = rank;
        process_log << "writing " << length << " integers" << endl;
        MPI_Barrier(writers);
        start = high_resolution_clock::now();
        IndividualProcessWrite((char*)(buffer.get()), length * sizeof(int32_t),
                               individual_file_path, process_log);
        MPI_Barrier(writers);
        end = high_resolution_clock::now();
        elapse = end - start;
        process_log << "individual processes writing outside of MPI:"
                    << elapse.count() << "ms" << endl;
        if (rank == 0)
            master_log << "writing individual file outside of MPI:"
                       << elapse.count() << "ms" << endl;
        boost::filesystem::remove(boost::filesystem::path(individual_file_path));
        string out_file_collective(out_file_base);
        out_file_collective.append(to_string(length)).append("_collective");
        if (mcp3d::IsFile(out_file_collective))
            remove(out_file_collective.c_str());
        // create file
        int error = MPI_File_open(
                writers, out_file_collective.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                info, &file);
        MPI_Barrier(writers);
        process_log << "file open" << endl;
        if (error)
            MPI_Abort(writers, mcp3d::MPI_FILE_CREATE_ERR);
        MPI_File_set_atomicity(file, false);
        // set file fiew
        MPI_Datatype arraytype;
        MPI_Type_contiguous(length, MPI_INT32_T, &arraytype);
        MPI_Type_commit(&arraytype);
        offset = long(writers_rank) * length * sizeof(int32_t);
        MPI_File_set_view(file, offset, MPI_INT32_T, arraytype, "native", info);
        MPI_Barrier(writers);
        process_log << "set view" << endl;
        MPI_Barrier(writers);
        // MPI write
        start = high_resolution_clock::now();
        //MPI_File_preallocate(
                //file, length * long(sizeof(int32_t)) * long(writers_size));
        //cout << "preallocate finish" << endl;
        MPI_File_write_all(file, buffer.get(), (int)length,
                           MPI_INT32_T, &status);
        cout << "write all finish" << endl;
        // MPI array type for file view
        end = high_resolution_clock::now();
        elapse = end - start;
        process_log << "file write:" << elapse.count() << "ms" << endl;
        MPI_Barrier(writers);
        end = high_resolution_clock::now();
        elapse = end - start;
        if (rank == 0)
            master_log << "writing " << length << " integers, time = "
                       << elapse.count() << " milliseconds\n" << endl;
        MPI_File_close(&file);
    }

    if (rank == 0)
        master_log.close();
    process_log.close();
    MPI_Finalize();
}

