//
// Created by muyezhu on 2/7/18.
//

#ifndef MCP3D_MPI_IO_UTIL_HPP
#define MCP3D_MPI_IO_UTIL_HPP

#include <fstream>
#include <mpi.h>
#include <unordered_map>
#include <nlohmann/json.hpp>

namespace mcp3d
{
enum ERROR_CODES
{
    MPI_FILE_CREATE_ERR
};

class RomioConfig
{
public:
    using json = nlohmann::json;
    explicit RomioConfig(const std::string& romio_config_path);
    const json& romio_config() { return romio_config_; }
    std::unordered_map<std::string, std::string> RomioHints()
    { return romio_config_.get<std::unordered_map<std::string, std::string>>(); };

private:
    std::string romio_json_path_;
    json romio_config_;
};

void SetMPIInfo(MPI_Info &info, const std::string& romio_json_path);

void LogMPIInfo(const MPI_Info &info, std::fstream &fs);

}
#endif //MCP3D_MPI_IO_UTIL_HPP
