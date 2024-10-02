//
// Created by muyezhu on 2/5/18.
//
#include <iostream>
#include <fstream>
#include <cstring>
#include <memory>
#include <cmath>
#include <chrono>
#include <hdf5.h>
#include <boost/filesystem.hpp>
#include "common/mcp3d_common.hpp"
#include "parallel/mpio_util.hpp"

using namespace std;

const int isilon_block_size = 8192;
const int default_chunk_dim = 1024;
const string trivial_dataset_name("trivial_data");
const string parallel_dataset_name("parallel_data");

string hdf5_parallel_test_out_path()
{
    return mcp3d::JoinPath({mcp3d::test_data_dir(), "infrastructure",
                            "hdf5_parallel"});
}

string trivial_hdf5_path()
{
    return mcp3d::JoinPath({hdf5_parallel_test_out_path(),
                            "trivial_create.hdf5"});
}

string parallel_hdf5_path()
{
    return mcp3d::JoinPath({hdf5_parallel_test_out_path(),
                            "parallel_create.hdf5"});
}

string romio_hints_path()
{
    return mcp3d::JoinPath({mcp3d::configs_dir(), "hdf5_romio_hints.json"});
}

void PrintFileChunkSize(long chunk_dim, long size, long bytes_per_elem)
{
    cout << "total processors = " << size << endl;
    cout << "chunk dimensions = [" << chunk_dim << ", " << chunk_dim << "]" << endl;
    cout << "bytes per element = " << bytes_per_elem << endl;
    long total_bytes = size * chunk_dim * chunk_dim * bytes_per_elem;
    double total_kb = total_bytes / 1024.0;
    double total_mb = total_kb / 1024.0;
    double total_gb = total_mb / 1024.0;
    cout << "file size = " << (total_gb > 1 ? to_string(total_gb) + "GB" : to_string(total_mb) + "MB") << endl;
}

bool TrivialCreate(uint8_t trivial_dataset_val)
{
    hid_t file, dataset, datatype, dataspace;
    herr_t status;
    int rank = 2;
    hsize_t n_elem;
    unique_ptr<hsize_t[]> dims(new hsize_t[rank]);
    dims[0] = 1024;
    dims[1] = 1024;
    n_elem =  dims[0] * dims[1];
    unique_ptr<uint8_t[]> data(new uint8_t[n_elem]);
    for (hsize_t i = 0; i < n_elem; ++i)
        data[i] = trivial_dataset_val;
    dataspace = H5Screate_simple(rank, dims.get(), nullptr);
    datatype = H5Tcopy(H5T_NATIVE_UCHAR);
    file = H5Fcreate(trivial_hdf5_path().c_str(), H5F_ACC_TRUNC,
                     H5P_DEFAULT, H5P_DEFAULT);
    if (file >= 0)
        cout << "created trial_create.hdf5" << endl;
    else MCP3D_TEST_ERROR("error creating trial_create.hdf5")
    dataset = H5Dcreate(file, trivial_dataset_name.c_str(), datatype, dataspace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    INIT_TIMER(0)
    TIC(0)
    status = H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data.get());
    TOC(0)
    REPORT_TIME_TO_COMPLETION("time to write dataset", 0)
    if (status >= 0)
        cout << "wrote 1024 * 1024 unsigned character to dataset" << endl;
    else MCP3D_TEST_ERROR("error writing to dataset")
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);
    int validate = H5Fis_hdf5(trivial_hdf5_path().c_str());
    if (validate < 0) MCP3D_TEST_ERROR("H5Fis_hdf5 function failure")
    else if (validate == 0) MCP3D_TEST_ERROR(
            trivial_hdf5_path() + " is not an hdf5 file")
    else
        cout << "write trivial hdf5 file success" << endl;
    return true;
}

bool TrivialRead(uint8_t trivial_dataset_val)
{
    hid_t file, dataset, dataspace, datatype, native_datatype;
    herr_t status;
    file = H5Fopen(trivial_hdf5_path().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) MCP3D_TEST_ERROR("file open failure")
    dataset = H5Dopen(file, trivial_dataset_name.c_str(), H5P_DEFAULT);
    if (dataset < 0) MCP3D_TEST_ERROR("dataset open failure")
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0) MCP3D_TEST_ERROR("dataspace open failure")
    int rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank < 0) MCP3D_TEST_ERROR("error reading rank of dataspace")
    datatype = H5Dget_type(dataset);
    native_datatype = H5Tget_native_type(datatype, H5T_DIR_ASCEND);
    if (!H5Tequal(native_datatype, H5T_NATIVE_UCHAR)) MCP3D_TEST_ERROR(
            "datatype incorrect")
    unique_ptr<hsize_t[]> dims(new hsize_t[rank]);
    rank = H5Sget_simple_extent_dims(dataspace, dims.get(), nullptr);
    if (rank < 0) MCP3D_TEST_ERROR("error reading dimensions of dataspace")
    hsize_t n_elem = 1;
    for (int i = 0; i < rank; ++i)
        n_elem *= dims[i];
    unique_ptr<uint8_t[]> data(new uint8_t[n_elem]);
    status = H5Dread(dataset, datatype, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, data.get());
    if (status < 0) MCP3D_TEST_ERROR("error reading dataset data")
    for (hsize_t i = 0; i < n_elem; ++i)
        if (data[i] != trivial_dataset_val) MCP3D_TEST_ERROR(
                "data value incorrect")
    cout << "read trivial hdf5 file success" << endl;
    H5Tclose(datatype);
    H5Tclose(native_datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);
    return true;
}

bool ParallelCreate(int chunk_dim)
{
    if (chunk_dim <= 0)
        chunk_dim = default_chunk_dim;
    const int block_dim = chunk_dim;
    int rank, size, is_mpi;
    MPI_Initialized(&is_mpi);
    if (!is_mpi) MCP3D_TEST_ERROR("MPI is not initialized")
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto dataset_val = (uint16_t)rank;
    auto block_xmax = (int)sqrt((double)size);
    auto block_ymax = block_xmax;
    int block_xid = rank % block_xmax;
    int block_yid = rank / block_xmax;

    // full array and chunk dimensions
    int n_dims = 2;
    hsize_t n_chunk_elem;
    unique_ptr<hsize_t[]> dims(new hsize_t[n_dims]);
    dims[0] = block_dim * hsize_t(block_ymax);
    dims[1] = block_dim * hsize_t(block_xmax);
    unique_ptr<hsize_t[]> chunk_dims(new hsize_t[n_dims]);
    chunk_dims[0] = (hsize_t)block_dim;
    chunk_dims[1] = (hsize_t)block_dim;

    if (rank == 0)
        PrintFileChunkSize(chunk_dim, size, 2);

    // define identifiers
    hid_t fcpl, fapl, file, dataset, datatype, file_space, mem_space, dcpl, dxpl;
    herr_t status;
    // prepare timers
    INIT_TIMER(0)
    // identical file creation property list
    TIC_WORLD(0)
    fcpl = H5Pcreate(H5P_FILE_CREATE);
    if (fcpl < 0) MCP3D_TEST_ERROR(
            "failed to create file creation property list")
    status = H5Pset_file_space_page_size(fcpl, isilon_block_size);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to file space page size property of file creation property list")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create file access property list", 0)

    // identical file access property list
    MPI_Info mpi_info;
    mcp3d::SetMPIInfo(mpi_info, romio_hints_path());
    TIC_WORLD(0)
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl < 0) MCP3D_TEST_ERROR("failed to create file access property list")
    status = H5Pset_driver(fapl, H5FD_MPIO, nullptr);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set driver property of file access property list")
    status = H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, mpi_info);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set mpio property of file access property list")
    status = H5Pset_coll_metadata_write(fapl, false);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set collective metadata write property of file access property list")
    status = H5Pset_alignment(fapl, 1, isilon_block_size);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set alignment property of file access property list")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create and set file access property list",0)

    // collective file creation
    TIC_WORLD(0)
    file = H5Fcreate(parallel_hdf5_path().c_str(), H5F_ACC_TRUNC,
                     H5P_DEFAULT, fapl);
    if (file < 0) MCP3D_TEST_ERROR("parallel hdf5 creation failure")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create file", 0)
    TIC_WORLD(0)
    status = H5Fset_mpi_atomicity(file, false);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set mpio atomicity for file")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "set file mpio atomicity", 0)

    // identical dataspace for all processors
    TIC_WORLD(0)
    file_space = H5Screate_simple(n_dims, dims.get(), nullptr);
    if (file_space < 0) MCP3D_TEST_ERROR("dataspace creation failure")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create dataspace", 0)

    // identical datatype for all processors
    TIC_WORLD(0)
    datatype = H5Tcopy(H5T_NATIVE_USHORT);
    if (datatype < 0) MCP3D_TEST_ERROR("datatype creation failure")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "datatype creation", 0)

    // identical dataset creation property list, allow chunks
    TIC_WORLD(0)
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0) MCP3D_TEST_ERROR(
            "failed to create dataset creation property list")
    status = H5Pset_layout(dcpl, H5D_CHUNKED);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set dataset layout to chunked for dataset creation property list")
    status = H5Pset_chunk(dcpl, n_dims, chunk_dims.get());
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set chunk size for dataset creation property list")
    status = H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_INCR);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set allocation time for dataset creation property list")
    status = H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set fill time for dataset creation property list")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create and set dataset creation property list",0)

    // collective dataset creation
    TIC_WORLD(0)
    dataset = H5Dcreate(file, parallel_dataset_name.c_str(), datatype, file_space,
                        H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (dataset < 0) MCP3D_TEST_ERROR("failed to create dataset")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create dataset", 0)

    // processor memory space
    TIC_WORLD(0)
    mem_space = H5Screate_simple(2, chunk_dims.get(), nullptr);
    // hyperslab set up
    unique_ptr<hsize_t[]> offset(new (nothrow) hsize_t[n_dims]);
    offset[0] = (hsize_t)block_yid * chunk_dims[0];
    offset[1] = (hsize_t)block_xid * chunk_dims[1];
    unique_ptr<hsize_t[]> count(new (nothrow) hsize_t[n_dims]);
    count[0] = 1;
    count[1] = 1;
    unique_ptr<hsize_t[]> stride(new (nothrow) hsize_t[n_dims]);
    stride[0] = 1;
    stride[1] = 1;
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
                                 offset.get(), stride.get(),
                                 count.get(), chunk_dims.get());
    if (status < 0) MCP3D_TEST_ERROR("failed to select hyperslab")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "select hyperslab", 0)

    // data buffer
    TIC_WORLD(0)
    n_chunk_elem =  chunk_dims[0] * chunk_dims[1];
    unique_ptr<uint16_t[]> data(new (nothrow)uint16_t[n_chunk_elem]);
    if (!data.get()) MCP3D_BAD_ALLOC("failed to allocate memory for data")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "allocate data buffer", 0)
    TIC_WORLD(0)
    uint16_t *data_ptr = data.get();
    for (hsize_t i = 0; i < n_chunk_elem; ++i)
        data_ptr[i] = dataset_val;
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "fill data buffer", 0)

    // data transfer property list
    TIC_WORLD(0)
    dxpl = H5Pcreate(H5P_DATASET_XFER);
    if (dxpl < 0) MCP3D_TEST_ERROR(
            "failed to create data transfer property list")
    status = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed to set mpio property for data transfer property list")
    status = H5Pset_dxpl_mpio_chunk_opt(dxpl, H5FD_MPIO_CHUNK_MULTI_IO);
    if (status < 0) MCP3D_TEST_ERROR(
            "failed set mpio chunk opt property for data transfer property list")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "create and set data transfer property list", 0)

    // write to file
    TIC_WORLD(0)
    status = H5Dwrite(dataset, H5T_NATIVE_USHORT,
                      mem_space, file_space, dxpl, data.get());
    if (status < 0) MCP3D_TEST_ERROR("failed to write to dataset")
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "write data", 0)

    // close handles
    TIC_WORLD(0)
    H5Dclose(dataset);
    H5Pclose(fapl);
    H5Pclose(fcpl);
    H5Pclose(dcpl);
    H5Pclose(dxpl);
    H5Sclose(file_space);
    H5Sclose(mem_space);
    H5Fclose(file);
    TOC_WORLD(0)
    MPI_REPORT_TIME_TO_COMPLETION(rank, "close handles", 0)

    double mem_usage = mcp3d::ProcPeakRAM("GB");
    if (rank == 0)
        cout << "peak memory usage: " << mem_usage << " GB" << endl;
    return true;
}

bool parallel_read(int chunk_dim)
{
    if (chunk_dim <= 0)
        chunk_dim = default_chunk_dim;
    const int block_dim = chunk_dim;
    int rank, size, is_mpi;
    MPI_Initialized(&is_mpi);
    if (!is_mpi) MCP3D_TEST_ERROR("MPI is not initialized")
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto dataset_val = (uint16_t)rank;
    auto block_xmax = (int)sqrt((double)size);
    auto block_ymax = block_xmax;
    int block_xid = rank % block_xmax;
    int block_yid = rank / block_xmax;


}

int main(int argc, char **argv)
{
    if (argc < 2 || (strcmp(argv[1], "trivial") != 0 && strcmp(argv[1], "parallel") != 0))
    {
        cout << "usage: ./hdf5_parallel_test trivial | parallel [chunk_dim]" << endl;
        exit(1);
    }
    using boost::filesystem::path;

    if (strcmp(argv[1], "trivial") == 0)
    {
        path pd(hdf5_parallel_test_out_path());
        path pf(trivial_hdf5_path());
        if (!boost::filesystem::is_directory(pd))
            boost::filesystem::create_directories(pd);
        if (boost::filesystem::exists(pf))
            boost::filesystem::remove(pf);
        uint8_t trivial_dataset_val = 255;
        TrivialCreate(trivial_dataset_val);
        TrivialRead(trivial_dataset_val);
    }
    else
    {
        MPI_Init(nullptr, nullptr);
        int size, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!IsInt(sqrt((double)size)) || size == 1)
            if (rank == 0)
            {
                cout << "please use number of processors equal to some square number greater than 1" << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        if (rank == 0)
        {
            path pd(hdf5_parallel_test_out_path());
            path pf(parallel_hdf5_path());
            if (!boost::filesystem::is_directory(pd))
                boost::filesystem::create_directories(pd);
            if (boost::filesystem::exists(pf))
                boost::filesystem::remove(pf);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int chunk_dim = -1;
        if (argc == 3)
        {
            try
            {
                chunk_dim = stoi(argv[2]);
            }
            catch(const invalid_argument &e)
            {
                cout << "usage: ./hdf5_parallel_test trivial | parallel [chunk_dim], "
                        "chunk_dim must be integer if given" << endl;
                exit(1);
            }
        }
        ParallelCreate(chunk_dim);
        parallel_read(chunk_dim);
        MPI_Finalize();
    }
}
