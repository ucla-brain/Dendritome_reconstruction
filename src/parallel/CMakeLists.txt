# mcp3d parallel module
set(MCP3D_PARALLEL_SRC
    mpio_util.hpp
    mpio_util.cpp)
add_library(mcp3d_parallel SHARED ${MCP3D_PARALLEL_SRC})

message("${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(write_image_pyramid_mpi write_image_pyramid_mpi.cpp)
set(WRITE_IMG_PYR_MPI_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/write_image_pyramid_mpi.sh
    ${CMAKE_CURRENT_SOURCE_DIR}/qsub_write_image_pyramid_mpi.sh)

add_executable(async_remove_path_mpi async_remove_path_mpi.cpp)
set(ASYNC_REMOVE_PATH_MPI_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/async_remove_path_mpi.sh
    ${CMAKE_CURRENT_SOURCE_DIR}/qsub_async_remove_path_mpi.sh)

add_executable(
   par_cc_3d
   par_connected_component_3d.hpp
   par_connected_component_3d.cpp
)
set(PAR_CC_3D_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/par_cc_3d.sh
    ${CMAKE_CURRENT_SOURCE_DIR}/qsub_par_cc_3d.sh)

target_link_libraries(mcp3d_parallel ${EXT_LINK_LIBS} mcp3d_common mcp3d_image mcp3d_algorithm)
target_link_libraries(write_image_pyramid_mpi ${EXT_LINK_LIBS} mcp3d_image)
target_link_libraries(async_remove_path_mpi ${EXT_LINK_LIBS} mcp3d_common)
target_link_libraries(par_cc_3d ${EXT_LINK_LIBS} mcp3d_common mcp3d_image mcp3d_algorithm)

add_custom_command(
        TARGET par_cc_3d POST_BUILD
        COMMAND cp ${PAR_CC_3D_FILES} ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(
        TARGET write_image_pyramid_mpi POST_BUILD
        COMMAND cp ${WRITE_IMG_PYR_MPI_FILES} ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(
        TARGET async_remove_path_mpi POST_BUILD
        COMMAND cp ${ASYNC_REMOVE_PATH_MPI_FILES} ${CMAKE_CURRENT_BINARY_DIR})