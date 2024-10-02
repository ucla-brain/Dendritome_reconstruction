//
// Created by muyezhu on 10/3/19.
//

#ifndef MCP3D_MCP3D_ALGORITHM_MACROS_HPP
#define MCP3D_MCP3D_ALGORITHM_MACROS_HPP

// note to not place brackets ({}) outside of using statements. these statements
// do not work properly when placed inside block scope
#define MCP3D_MGRAPH_ALIASES                                    \
    using typename mcp3d::MGraph<Vertex>::EdgeType;             \
    using mcp3d::MGraph<Vertex>::StoreEdges;                    \
    using mcp3d::MGraph<Vertex>::Size;                          \
    using mcp3d::MGraph<Vertex>::HasVertex;                     \
    using mcp3d::MGraph<Vertex>::At;                            \
    using mcp3d::MGraph<Vertex>::operator[];                    \
    using mcp3d::MGraph<Vertex>::Vertices;                      \

#endif //MCP3D_MCP3D_ALGORITHM_MACROS_HPP
