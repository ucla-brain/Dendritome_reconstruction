set(MCP3D_ALGORITHM_SRC
    mcp3d_algorithm_macros.hpp
    mcp3d_connected_component_3d.hpp
    mcp3d_vertex.hpp
    mcp3d_graph.hpp
    mcp3d_heap.hpp
    mcp3d_dijkstra.hpp
    mcp3d_all_path_pruning.hpp
    mcp3d_connected_component_3d.cpp)

add_library(mcp3d_algorithm SHARED ${MCP3D_ALGORITHM_SRC})
target_link_libraries(mcp3d_algorithm ${EXT_LINK_LIBS} mcp3d_common mcp3d_image)

add_executable(mcp3d_all_path_pruning  ${MCP3D_ALGORITHM_SRC} mcp3d_all_path_pruning.cpp)
target_link_libraries(mcp3d_all_path_pruning mcp3d_algorithm)

install(TARGETS mcp3d_algorithm LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
set_target_properties(mcp3d_algorithm PROPERTIES INSTALL_RPATH \$ORIGIN)

install(TARGETS mcp3d_all_path_pruning RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
set_target_properties(mcp3d_all_path_pruning PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)