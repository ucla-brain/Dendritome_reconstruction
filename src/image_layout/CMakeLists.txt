set(MCP3D_IMAGE_LAYOUT_SRC
    mcp3d_image_layout.hpp
    mcp3d_image_layout_constants.hpp
    mcp3d_channel_layout.hpp
    mcp3d_channel_pyr_slices.hpp
    mcp3d_volume_layout.hpp
    mcp3d_volume_metadata.hpp
    mcp3d_channel_layout.cpp
    mcp3d_channel_pyr_slices.cpp
    mcp3d_volume_layout.cpp
    mcp3d_volume_metadata.cpp)

add_library(mcp3d_image_layout SHARED ${MCP3D_IMAGE_LAYOUT_SRC})
target_link_libraries(mcp3d_image_layout
                      ${EXT_LINK_LIBS}
                      mcp3d_common
                      mcp3d_image_interface)

install(TARGETS mcp3d_image_layout LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
set_target_properties(mcp3d_image_layout PROPERTIES INSTALL_RPATH \$ORIGIN)