//
// Created by mzhu on 1/24/18.
//
#include <boost/predef.h>
#include <boost/filesystem.hpp>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_paths.hpp"

using namespace std;
using namespace boost;

char mcp3d::SeparatorChar()
{
    #if BOOST_OS_WINDOWS
        return '\\';
    #else
        return '/';
    #endif
}

char mcp3d::SeparatorChar(bool is_linux)
{
    if (is_linux)
        return '/';
    else
        return '\\';
}

string mcp3d::src_dir()
{
    filesystem::path self_path(__FILE__);
    string d(self_path.parent_path().parent_path().string());
    if (!filesystem::is_directory(filesystem::path(d)))
        MCP3D_RUNTIME_ERROR(d + " does not exist")
    return d;
}

string mcp3d::parallel_module_dir()
{
    string d(mcp3d::src_dir());
    d.append(1, mcp3d::SeparatorChar());
    d.append("parallel");
    if (!filesystem::is_directory(filesystem::path(d)))
        MCP3D_RUNTIME_ERROR(d + " does not exist")
    return d;
}

string mcp3d::test_data_dir()
{
    filesystem::path src(mcp3d::src_dir());
    string d(src.parent_path().string());
    d.append(1, mcp3d::SeparatorChar());
    d.append("test_data");
    d.append(1, mcp3d::SeparatorChar());
    d.append("src");
    if (!filesystem::is_directory(filesystem::path(d)))
        MCP3D_RUNTIME_ERROR(d + " does not exist")
    return d;
}

string mcp3d::benchmark_module_dir()
{
    string d(mcp3d::src_dir());
    d.append(1, mcp3d::SeparatorChar());
    d.append("benchmark");
    if (!filesystem::is_directory(filesystem::path(d)))
        MCP3D_RUNTIME_ERROR(d + " does not exist")
    return d;
}
