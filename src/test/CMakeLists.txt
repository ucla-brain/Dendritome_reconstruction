# test module general
set(UNIT_TEST_SRC
    test_cc_label_3d_single_core.hpp
    test_cc_label_3d_single_core.cpp
    test_mcp3d_exceptions.cpp
    test_utility_internal.hpp
    test_utility_internal.cpp
    test_mcp3d_utility.cpp
    test_openmp.cpp
    test_mcp3d_voxel_types.cpp
    test_mcp3d_channel_pyr_slices.hpp
    test_mcp3d_channel_pyr_slices.cpp
    test_mcp3d_channel_layout.cpp


    #test_mcp3d_image_utils.cpp
    #test_mcp3d_image_maths.cpp
    #test_mcp3d_image_info.cpp
    #test_mcp3d_image_view.cpp
    #test_mcp3d_image.cpp
    #test_mcp3d_tiff_format.cpp
    #test_mcp3d_ome_tiff_format.cpp
    #test_mcp3d_hdf5_format.cpp
    #test_mcp3d_vertex.cpp
    #test_mcp3d_heap.cpp
    #test_mcp3d_graph.cpp
    #test_mcp3d_dijkstra.hpp
    #test_mcp3d_dijkstra.cpp
    #test_mcp3d_all_path_pruning.cpp
    mcp3d_unit_tests.cpp)
add_executable(mcp3d_unit_tests ${UNIT_TEST_SRC})
target_link_libraries(mcp3d_unit_tests ${EXT_LINK_LIBS} mcp3d_common mcp3d_image_layout mcp3d_image mcp3d_algorithm)
set_target_properties(mcp3d_unit_tests PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)

if(MPI)
    set(MPI_TEST_SRC
        test_mpi.cpp
        test_utility_internal.hpp
        test_utility_internal.cpp
        test_async_remove_path_mpi.cpp
        mcp3d_mpi_tests.cpp)
    add_executable(mcp3d_mpi_tests ${MPI_TEST_SRC})
    target_link_libraries(mcp3d_mpi_tests ${EXT_LINK_LIBS} mcp3d_common mcp3d_image_layout mcp3d_image mcp3d_algorithm)
    set_target_properties(mcp3d_mpi_tests PROPERTIES INSTALL_RPATH \$ORIGIN/../lib)
endif()


install(TARGETS mcp3d_unit_tests
        RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
if(MPI)
    install(TARGETS mcp3d_mpi_tests
            RUNTIME DESTINATION ${MCP3D_INSTALL_DIR}/bin)
endif()