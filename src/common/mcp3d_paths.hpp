//
// Created by mzhu on 1/24/18.
//

#ifndef MCP3D_MCP3D_PATHS_HPP
#define MCP3D_MCP3D_PATHS_HPP

namespace mcp3d
{
    char SeparatorChar();

    char SeparatorChar(bool is_linux);

    inline std::string mcp_dir() { return "/ifs/loni/faculty/dong/mcp"; }

    inline std::string z_drive()  { return "Z:\\"; }

    std::string src_dir();

    std::string parallel_module_dir();

    std::string test_data_dir();

    std::string benchmark_module_dir();
}

#endif //MCP3D_MCP3D_PATHS_HPP
