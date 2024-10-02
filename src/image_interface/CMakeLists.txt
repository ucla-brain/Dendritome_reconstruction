set(MCP3D_IMAGE_INTERFACE_SRC
    mcp3d_voxel_types.hpp
    mcp3d_file_formats.hpp
    mcp3d_tiff_utils.hpp
    mcp3d_tiff_io.hpp
    mcp3d_hdf5_utils.hpp
    mcp3d_imaris_util.hpp
    mcp3d_voxel_types.cpp
    mcp3d_file_formats.cpp
    mcp3d_tiff_utils.cpp
    mcp3d_tiff_io.cpp
    mcp3d_hdf5_utils.cpp
    mcp3d_imaris_util.cpp)

add_library(mcp3d_image_interface SHARED ${MCP3D_IMAGE_INTERFACE_SRC})
target_link_libraries(mcp3d_image_interface
                      ${EXT_LINK_LIBS}
                      mcp3d_common)

install(TARGETS mcp3d_image_interface LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
set_target_properties(mcp3d_image_interface PROPERTIES INSTALL_RPATH \$ORIGIN)
