#!/bin/bash

script_dir="$(dirname "${BASH_SOURCE[0]}")"
project_dir="$(bash "$script_dir/get_project_dir.sh")"
cd "$project_dir"
current_branch=$(git rev-parse --abbrev-ref HEAD)
echo "project directory: $project_dir"
build_dir="$project_dir/grid_build_""$current_branch"
if [[ "$1" = "clean-build" ]]
then
    echo "removing $build_dir"
    rm -rf "$build_dir"
fi
echo "build directory: $build_dir"
cmake_bin="/ifs/loni/faculty/dong/mcp/utils/cmake3.8/cmake-3.8.2-Linux-x86_64/bin/cmake"
# you need to have PATH and LD_LIBRARY variabies with gcc-7.2 and mpich-3.1.4
# entries prepending Centos system default
export CC=/usr/local/gcc-7.2.0/bin/gcc
export CXX=/usr/local/gcc-7.2.0/bin/g++
export MPI_C=/usr/local/mpich-3.1.3/bin/mpicc
export MPI_CXX=/usr/local/mpich-3.1.3/bin/mpicxx
export MPIEXEC=/usr/local/mpich-3.1.3/bin/mpiexec

mkdir -p "$build_dir"
cd "$build_dir"
"$cmake_bin" -DCMAKE_BUILD_TYPE=Release -DMPI=TRUE "$project_dir"'/src'
make
exit 0
