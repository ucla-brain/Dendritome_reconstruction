//
// Created by muyezhu on 9/15/17.
//

#ifndef MCP3D_MCP3D_TYPES_HPP
#define MCP3D_MCP3D_TYPES_HPP

#include <cstdint>
#include <type_traits>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

namespace mcp3d
{
template<typename T>
using MCPTensor0D = Eigen::Tensor<T, 0, Eigen::RowMajor + Eigen::AutoAlign>;

template<typename T>
using MCPTensor1D = Eigen::Tensor<T, 1, Eigen::RowMajor + Eigen::AutoAlign>;

template<typename T>
using MCPTensor1DMap = Eigen::TensorMap<Eigen::Tensor<T, 1, Eigen::RowMajor>, Eigen::Aligned>;

template<typename T>
using MCPTensor2D = Eigen::Tensor<T, 2, Eigen::RowMajor + Eigen::AutoAlign>;

template<typename T>
using MCPTensor2DMap = Eigen::TensorMap<Eigen::Tensor<T, 2, Eigen::RowMajor>, Eigen::Aligned>;

template<typename T>
using MCPTensor3D = Eigen::Tensor<T, 3, Eigen::RowMajor + Eigen::AutoAlign>;

template<typename T>
using MCPTensor3DMap = Eigen::TensorMap<Eigen::Tensor<T, 3, Eigen::RowMajor>, Eigen::Aligned>;



}
#endif //MCP3D_MCP3D_TYPES_HPP
