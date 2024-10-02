message("external tool: spatial distance (from Vaa3D)")
set(VAA3D_SPATIAL_DISTANCE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/../xyz.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../neuron_tree.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../swc_file.hpp
    v_neuronswc.hpp
    neuron_sim_scores.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../swc_file.cpp
    v_neuronswc.cpp
    neuron_sim_scores.cpp
    spatial_distance.cpp
    )

add_executable(spatial_distance ${VAA3D_SPATIAL_DISTANCE_SRC})
target_link_libraries(spatial_distance ${EXT_LINK_LIBS} mcp3d_common)

install(TARGETS spatial_distance RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
set_target_properties(spatial_distance PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)