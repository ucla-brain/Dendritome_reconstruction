//
// Created by muyezhu on 2/7/18.
//
#include <unordered_map>
#include "boost/filesystem.hpp"
#include "common/mcp3d_common.hpp"
#include "parallel/mpio_util.hpp"

using namespace std;

void mcp3d::SetMPIInfo(MPI_Info &info, const string& romio_json_path)
{
    RomioConfig romio_config(romio_json_path);
    unordered_map<string, string> romio_hints = romio_config.RomioHints();
    MPI_Info_create(&info);
    for (auto hint: romio_hints)
        MPI_Info_set(info, (hint.first).c_str(), (hint.second).c_str());
}

void mcp3d::LogMPIInfo(const MPI_Info &info, fstream &fs)
{
    int n_keys, value_len, flag;
    char key[255];
    char value[255];
    MPI_Info_get_nkeys(info, &n_keys);
    fs << "MPI_Info: " << endl;
    for (int i = 0; i < n_keys; ++i)
    {
        MPI_Info_get_nthkey(info, i, key);
        MPI_Info_get_valuelen(info, key, &value_len, &flag);
        // MPICH and OpenMPI has 1 character difference in value_len
        // use + 1 on cluster.
        MPI_Info_get(info, key, value_len + 1, value, &flag);
        if (flag)
            fs << "key:" << key << ", " << "value:" << string(value) << endl;
        else
            fs << i << "th key can not be retrieved" << endl;
    }
}


mcp3d::RomioConfig::RomioConfig(const string& romio_config_path):
                                romio_json_path_(romio_config_path),
                                romio_config_(json {})
{
    using namespace boost;
    if (!filesystem::is_regular_file(filesystem::path(romio_json_path_)))
    MCP3D_INVALID_ARGUMENT(romio_json_path_ + " does not exist")
    ifstream ifs(romio_json_path_);
    ifs >> romio_config_;
    ifs.close();
}