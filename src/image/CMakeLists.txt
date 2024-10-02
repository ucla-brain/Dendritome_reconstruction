set(MCP3D_IMAGE_SRC
    mcp3d_image_macros.hpp
    mcp3d_image_constants.hpp
    mcp3d_image_utils.hpp
    mcp3d_image_maths.hpp
    mcp3d_channel_info.hpp
    mcp3d_image_info.hpp
    mcp3d_image_view.hpp
    mcp3d_image_formats.hpp
    mcp3d_tiff_format.hpp
    mcp3d_ome_tiff_format.hpp
    mcp3d_imaris_format.hpp
    mcp3d_hdf5_format.hpp
    mcp3d_image_io.hpp
    mcp3d_image_base.hpp
    mcp3d_image.hpp
    mcp3d_image_in_memory.hpp
    mcp3d_image_voxel_iterator.hpp
    mcp3d_image_utils.cpp
    mcp3d_channel_info.cpp
    mcp3d_image_info.cpp
    mcp3d_image_view.cpp
    mcp3d_image_formats.cpp
    mcp3d_tiff_format.cpp
    mcp3d_ome_tiff_format.cpp
    mcp3d_imaris_format.cpp
    mcp3d_hdf5_format.cpp
    mcp3d_image_io.cpp
    mcp3d_image_base.cpp
    mcp3d_image.cpp
    mcp3d_image_in_memory.cpp)

add_library(mcp3d_image SHARED ${MCP3D_IMAGE_SRC})
target_link_libraries(mcp3d_image
                      ${EXT_LINK_LIBS}
                      mcp3d_common
                      mcp3d_image_interface
                      mcp3d_image_layout)

install(TARGETS mcp3d_image LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
set_target_properties(mcp3d_image PROPERTIES INSTALL_RPATH \$ORIGIN)
