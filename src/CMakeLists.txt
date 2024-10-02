cmake_minimum_required(VERSION 3.5)
project(mcp3d)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ENV{OMP_CANCELLATION} true)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message("cmake build type: ${CMAKE_BUILD_TYPE}")

# -ffast-math compiles slower code
# use readelf -p .GCC.command.line executable to see compiler flags
# -fprofile-arcs -ftest-coverage options for gcov to generate code coverage information
# -pg: flat for gprof. execute the binary will produce gmon.out. then use grof program_name -> out_file to get flat profile and call graph
# only add -pg if profiling needed
if (PROFILE)
    message("compiling with -pg flag for gprof. python binding will not work")
    set(CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS} -fopenmp -pthread -frecord-gcc-switches -fprofile-arcs -ftest-coverage -march=native -Wall -g -pg")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fopenmp -pthread -frecord-gcc-switches -fprofile-arcs -ftest-coverage -std=c++14 -march=native -Wall -g -pg")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -fopenmp -pthread -frecord-gcc-switches -march=native  -Wall -g -pg -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -fopenmp -pthread -frecord-gcc-switches -std=c++11 -march=native -Wall -g -pg -O3  -DNDEBUG")
else()
    set(CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS} -fopenmp -pthread -frecord-gcc-switches -fprofile-arcs -ftest-coverage -march=native -Wall -g")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fopenmp -pthread -frecord-gcc-switches -fprofile-arcs -ftest-coverage -std=c++14 -march=native -Wall -g")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -fopenmp -pthread -frecord-gcc-switches -march=native -fomit-frame-pointer -Wall -g -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -fopenmp -pthread -frecord-gcc-switches -std=c++11 -march=native -Wall -g -O3 -fomit-frame-pointer -DNDEBUG")
endif()

if (VERBOSE)
    add_definitions(-DVERBOSE)
    message("verbose mode")
endif()

# host dependent settings
cmake_host_system_information(RESULT MACHINE QUERY HOSTNAME)
if (${MACHINE} STREQUAL "c2001")
    # add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
    # comment out openblas build for now. 20% slower than eigen in convolution
    message("building on ${MACHINE}...")
    set(MPI true)
    set(CMAKE_CC_COMPILER "/usr/local/gcc-7.2.0/bin/gcc")
    set(CMAKE_CXX_COMPILER "/usr/local/gcc-7.2.0/bin/g++")
    set(MPI_C_COMPILER "/usr/local/mpich-3.1.3/bin/mpicc")
    set(MPI_CXX_COMPILER "/usr/local/mpich-3.1.3/bin/mpicxx")
    set(MPIEXEC "/usr/local/mpich-3.1.3/bin/mpiexec")
    include_directories("/ifs/loni/faculty/dong/mcp/utils/libtiff-4.0.8/include")
    link_directories("/ifs/loni/faculty/dong/mcp/utils/libtiff-4.0.8/lib")
else()
    message("building on local machine...")
endif()

include_directories(${PROJECT_SOURCE_DIR})

if(MPI)
    add_definitions(-DSUPPORT_MPI)
    find_package(MPI REQUIRED)
    message("building with MPI support...")
    message(STATUS "MPI_INCLUDE_DIRS: ${MPI_INCLUDE_PATH}")
    message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
else()
    message(STATUS "not building MPI components")
endif()

set(3RD_PARTY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party)

find_package(EIGEN REQUIRED)
message(STATUS "EIGEN_INCLUDE_PATH: ${EIGEN_INCLUDE_PATH}")

find_package(GTEST REQUIRED)
message(STATUS "GTEST_INCLUDE_PATH: ${GTEST_INCLUDE_PATH}")
message(STATUS "GTEST_LIBRARIES: ${GTEST_LIBRARIES}")

find_package(NLOHMANN REQUIRED
             NAMES nlohmann_json
             PATHS ${3RD_PARTY_PATH}/nlohmann
             NO_DEFAULT_PATH)
message(
STATUS
"NLOHMANN package include directory generated with expression \
$<TARGET_PROPERTY:nlohmann_json,INTERFACE_INCLUDE_DIRECTORIES>")

find_package(BOOST REQUIRED)
message(STATUS "BOOST_INCLUDE_PATH: ${BOOST_INCLUDE_PATH}")
message(STATUS "BOOST_LIBRARY_PATH: ${BOOST_LIBRARY_PATH}")

find_package(HDF5 COMPONENTS C REQUIRED shared
             PATHS ${3RD_PARTY_PATH}/hdf5/share
             NO_DEFAULT_PATH)
message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")
message(STATUS "HDF5 TARGETS: ${HDF5_C_SHARED_LIBRARY}")
if (${HDF5_ENABLE_PARALLEL})
    message("warning: parallel hdf5 library used. non MPI build will not suceed")
endif()

find_package(OpenCV REQUIRED core imgcodecs imgproc
             PATHS ${3RD_PARTY_PATH}/opencv
             NO_DEFAULT_PATH)
message(STATUS "OPENCV_INCLUDE_DIR: ${OpenCV_INCLUDE_DIRS}")
message(STATUS "OPENCV_LIBRARYS: ${OpenCV_LIBS}")

set(EXT_INCLUDE_PATHS
    ${MPI_INCLUDE_PATH}
    ${EIGEN_INCLUDE_PATH}
    ${GTEST_INCLUDE_PATH}
    # see nlohmann_jsonTargets.cmake target property @line 50
    $<TARGET_PROPERTY:nlohmann_json,INTERFACE_INCLUDE_DIRECTORIES>
    ${BOOST_INCLUDE_PATH}
    ${HDF5_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS})
include_directories(${EXT_INCLUDE_PATHS})

link_directories(${BOOST_LIBRARY_PATH})
set(BOOST_LIBRARIES
    ${BOOST_LIBRARY_PATH}/libboost_serialization.so
    ${BOOST_LIBRARY_PATH}/libboost_filesystem.so
    ${BOOST_LIBRARY_PATH}/libboost_system.so)

set(HDF5_LIBRARIES
    ${HDF5_C_SHARED_LIBRARY})

set(EXT_LINK_LIBS
    libm.so
    libtiff.so
    libpthread.so
    liblapack.so
    ${MPI_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${BOOST_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${OpenCV_LIBS})
message("external link libraries: ${EXT_LINK_LIBS}")

# eigen gives ignored attributes warnings in compilation
# http://eigen.tuxfamily.org/bz_attachmentbase/attachment.cgi?id=708
# seems harmless

# install boost libraries
get_filename_component(MCP3D_INSTALL_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
foreach(BOOST_LIBRARY_LINK ${BOOST_LIBRARIES})
    get_filename_component(BOOST_LIBRARY ${BOOST_LIBRARY_LINK} REALPATH)
    install(FILES ${BOOST_LIBRARY}
            DESTINATION ${MCP3D_INSTALL_DIR}/lib)
endforeach()

# install hdf5 libraries. find soname and realpath of target, install both
# see 3rd_party/hdf5/share/cmake/hdf5-targets-release.cmake:21
set(HDF5_LIBRARY_PATH ${3RD_PARTY_PATH}/hdf5/lib)
get_target_property(HDF5_C_SHARED_LIBRARY_SONAME
                    ${HDF5_C_SHARED_LIBRARY} IMPORTED_SONAME_RELEASE)
install(FILES ${HDF5_LIBRARY_PATH}/${HDF5_C_SHARED_LIBRARY_SONAME}
        DESTINATION ${MCP3D_INSTALL_DIR}/lib)
get_filename_component(HDF5_C_SHARED_LIBRARY_LOCATION
                      ${HDF5_LIBRARY_PATH}/${HDF5_C_SHARED_LIBRARY_SONAME}
                      REALPATH)
install(FILES ${HDF5_C_SHARED_LIBRARY_LOCATION}
        DESTINATION ${MCP3D_INSTALL_DIR}/lib)

# install opencv libraries
foreach(OpenCV_LIB ${OpenCV_LIBS})
    get_target_property(OpenCV_LIB_REAL_PATH ${OpenCV_LIB} IMPORTED_LOCATION_RELEASE)
    install(FILES ${OpenCV_LIB_REAL_PATH} DESTINATION ${MCP3D_INSTALL_DIR}/lib)
    get_filename_component(OpenCV_LIB_DIR ${OpenCV_LIB_REAL_PATH} DIRECTORY)
    get_target_property(OpenCV_LIB_SONAME ${OpenCV_LIB} IMPORTED_SONAME_RELEASE)
    install(FILES ${OpenCV_LIB_DIR}/${OpenCV_LIB_SONAME} DESTINATION ${MCP3D_INSTALL_DIR}/lib)
endforeach()

# project target install will set target INSTALL_RPATH property with \$ORIGIN
# to produce relocatable libraries and executables
add_subdirectory(common)
add_subdirectory(image_interface)
add_subdirectory(image_layout)
add_subdirectory(image)
add_subdirectory(algorithm)
add_subdirectory(test)
add_subdirectory(external_tools)
# python wrapper
add_subdirectory(python)

if(MPI)
    add_subdirectory(parallel)
    if (PROFILE)
        add_subdirectory(benchmark)
    endif()
    add_subdirectory(test/infrastructure)
endif()
