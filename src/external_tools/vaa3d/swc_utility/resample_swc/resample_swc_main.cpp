/* resample_swc_func.cpp
 * This is a plugin to resample neuron swc subject to a fixed step length.
 * 2012-03-02 : by Yinan Wan
 */

#include <cstdint>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_utility.hpp"
#include "vaa3d/swc_utility/swc_file.hpp"
#include "resample_swc.hpp"

using namespace std;

void PrintUsage()
{
    cout << "(version 1.0) Resample points in a swc file subject to a fixed step length. Developed by Yinan Wan 2012-03-02" << endl;
    cout << "usage: ./resample_swc <input_swc_path> [step] [output_swc_path]" << endl;
    cout << "<input_swc_path>        input file path" << endl;
    cout << "[step_length]           step length for resampling, default to 10.0" << endl;
    cout << "[output_swc_path]       output file path, default to <input_swc_basename>_resampled.swc" << endl;
}

int main(int argc, char** argv)
{
	if (argc < 2)
    {
        PrintUsage();
        MCP3D_MESSAGE("too few arguments. exiting")
        return 1;
    }

	string input_swc_path(argv[1]);
    if (!mcp3d::IsFile(input_swc_path))
    {
        MCP3D_MESSAGE(input_swc_path + " is not a file")
        return 1;
    }

    double step = 10.0;
    if (argc >= 3)
    {
        try
        {
            step = mcp3d::ParseDouble(argv[2]);
        }
        catch (const exception& e)
        {
            MCP3D_MESSAGE(e.what())
            return 1;
        }
    }

    string output_swc_path = argc >= 4 ?
                             argv[3] : input_swc_path.substr(0, input_swc_path.size() - 4) + "_resample.swc";

	mcp3d::ResampleSwc(input_swc_path, step, output_swc_path);

	return 0;
}
