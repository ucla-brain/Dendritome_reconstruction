message("external tool: sort swc (from Vaa3D)")
set(SORT_SWC_SRC sort_swc.h sort_swc_main.cpp)

add_executable(sort_swc ${SORT_SWC_SRC})
target_link_libraries(sort_swc ${EXT_LINK_LIBS} mcp3d_common)

install(TARGETS sort_swc RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
set_target_properties(sort_swc PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)

