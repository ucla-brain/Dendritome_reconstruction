//resample neuronTree subject to a step length
//2012-02-29 by Yinan Wan
//2012-03-05 Yinan Wan: interpolate radius
#ifndef MCP3D_VAA3D_RESAMPLING_HPP
#define MCP3D_VAA3D_RESAMPLING_HPP

#include <vector>
#include <cmath>
#include "vaa3d/swc_utility/neuron_tree.hpp"

NeuronTree ResampleNeuronTree(NeuronTree input, double step);

namespace mcp3d
{
void ResampleSwc(const std::string& input_swc_path, double step = 10.0, const std::string& output_swc_path = std::string());
}


#endif  // MCP3D_VAA3D_RESAMPLING_HPP
