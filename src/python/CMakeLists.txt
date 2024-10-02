if(PYTHON)
    # temperarily hard code Python.h.
    if(NOT PYTHON_VERSION)
        set(PYTHON_VERSION "2")
    endif()
    if("${PYTHON_VERSION}" STREQUAL "2")
        find_package(PythonInterp 2.7 REQUIRED)
        find_package(PythonLibs 2.7 REQUIRED)
    # if no minor python 3 version given, use the highest available
    elseif("${PYTHON_VERSION}" STREQUAL "3")
        find_package(PythonInterp 3.7)
        if(NOT PYTHONINTERP_FOUND)
            find_package(PythonInterp 3.6)
            if(NOT PYTHONINTERP_FOUND)
                find_package(PythonInterp 3.5)
            endif()
        endif()
    elseif("${PYTHON_VERSION}" STREQUAL "3.7")
        find_package(PythonInterp 3.7 REQUIRED)
        find_package(PythonLibs 3.7 REQUIRED)
    elseif("${PYTHON_VERSION}" STREQUAL "3.6")
        find_package(PythonInterp 3.6 REQUIRED)
        find_package(PythonLibs 3.6 REQUIRED)
    elseif("${PYTHON_VERSION}" STREQUAL "3.5")
        find_package(PythonInterp 3.5 REQUIRED)
        find_package(PythonLibs 3.5 REQUIRED)
    endif()
    if(NOT PYTHONINTERP_FOUND)
        message(FATAL_ERROR "python interpreter not found")
    endif()
    find_package(PythonLibs REQUIRED)

    message("python version: ${PYTHON_VERSION_STRING}")
    message("PYTHON_INCLUDE_DIR: ${PYTHON_INCLUDE_PATH}")
    message("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
    message("cython version known to work: 0.29.6")


    set(PYTHON_CMAKE_CXX_FLAGS "-fopenmp -pthread -frecord-gcc-switches -fno-strict-aliasing -std=c++11 -march=native -Wall -g -O2")
    set(CYTHON_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/mcp3d_cython.pxd
        ${CMAKE_CURRENT_SOURCE_DIR}/pymcp3d.pyx)
    if ("${PYTHON_VERSION}" STREQUAL "2")

        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d.cpp
                           COMMAND cython -a -2 --cplus ${CYTHON_SRC} -o ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d.cpp
                           DEPENDS ${CYTHON_SRC})
        add_library(pymcp3d SHARED ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d.cpp)
        # -fno-strict-aliasing to suppress warnings
        # https://github.com/cython/cython/issues/958 https://github.com/boostorg/python/issues/179
        set_target_properties(pymcp3d PROPERTIES PREFIX "" COMPILE_FLAGS "${PYTHON_CMAKE_CXX_FLAGS}")
        target_include_directories(pymcp3d PUBLIC ${PYTHON_INCLUDE_PATH})
        target_link_libraries(pymcp3d ${PYTHON_LIBRARIES} mcp3d_common mcp3d_image mcp3d_algorithm)
    else()
        set(CYTHON3_SRC
            ${CMAKE_CURRENT_SOURCE_DIR}/mcp3d_cython.pxd
            ${CMAKE_CURRENT_SOURCE_DIR}/pymcp3d3.pyx)
        # copy pymcp3d.pyx to pymcp3d3.pyx. cython needs the pyx and cpp files
        # to have the same name for current PyInit_ exports
        add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/pymcp3d3.pyx
                                  ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d3.cpp
                           COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/pymcp3d.pyx ${CMAKE_CURRENT_SOURCE_DIR}/pymcp3d3.pyx
                           COMMAND cython -a -3 --cplus ${CYTHON3_SRC} -o ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d3.cpp
                           DEPENDS ${CYTHON_SRC})
        add_library(pymcp3d3 SHARED ${CMAKE_CURRENT_BINARY_DIR}/pymcp3d3.cpp)
        target_include_directories(pymcp3d3 PUBLIC ${PYTHON_INCLUDE_PATH})
        target_link_libraries(pymcp3d3 ${PYTHON_LIBRARIES} mcp3d_common mcp3d_image mcp3d_algorithm)
        set_target_properties(pymcp3d3 PROPERTIES PREFIX "" COMPILE_FLAGS "${PYTHON_CMAKE_CXX_FLAGS}")
    endif()

    if ("${PYTHON_VERSION}" STREQUAL "2")
        install(TARGETS pymcp3d LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
        set_target_properties(pymcp3d PROPERTIES INSTALL_RPATH \$ORIGIN)
    else()
        install(TARGETS pymcp3d3 LIBRARY DESTINATION ${MCP3D_INSTALL_DIR}/lib)
        set_target_properties(pymcp3d3 PROPERTIES INSTALL_RPATH \$ORIGIN)
    endif()
else()
    message(STATUS "not compiling python bindings")
endif()