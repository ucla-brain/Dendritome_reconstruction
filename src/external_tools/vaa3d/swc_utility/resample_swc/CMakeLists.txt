message("external tool: resample swc (from Vaa3D)")
set(RESAMPLE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/../swc_file.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../neuron_tree.hpp
    resample_swc.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../swc_file.cpp
    resample_swc.cpp
    resample_swc_main.cpp)

add_executable(resample_swc ${RESAMPLE_SRC})
target_link_libraries(resample_swc ${EXT_LINK_LIBS} mcp3d_common)

install(TARGETS resample_swc RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
set_target_properties(resample_swc PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)