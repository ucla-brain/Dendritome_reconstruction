message("external tool: all path pruning2 (from Vaa3D)")
set(VAA3D_APP2_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/markers.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/swc_file.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/resample_swc/resample_swc.hpp
    fastmarching_dt.h
    fastmarching_macro.h
    fastmarching_tree.h
    fastmarching_general.hpp
    heap.h
    hierarchy_prune.h
    marker_radius.h
    smooth_curve.h
    upwind_solver.h
    app2_parameters.hpp
    app2_image_preprocess.hpp
    all_path_pruning2.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/markers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/swc_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../../swc_utility/resample_swc/resample_swc.cpp
    app2_parameters.cpp
    app2_image_preprocess.cpp
    all_path_pruning2.cpp
    app2_main.cpp)
add_executable(vaa3d_app2 ${VAA3D_APP2_SRC})
target_link_libraries(vaa3d_app2 ${EXT_LINK_LIBS} mcp3d_common mcp3d_image)

install(TARGETS vaa3d_app2 RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
set_target_properties(vaa3d_app2 PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)