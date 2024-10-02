# test module for cluster infrastructure
add_executable(mpi_simple_test mpi_simple_test.cpp)
add_executable(mpi_io_test mpi_io_test.cpp)

target_link_libraries(mpi_simple_test ${EXT_LINK_LIBS} mcp3d_common mcp3d_image mcp3d_algorithm mcp3d_parallel)
target_link_libraries(mpi_io_test ${EXT_LINK_LIBS} mcp3d_common mcp3d_image mcp3d_algorithm mcp3d_parallel)