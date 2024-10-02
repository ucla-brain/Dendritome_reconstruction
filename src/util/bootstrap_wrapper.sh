#!/usr/bin/env bash

script_dir="$(dirname "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "get_project_dir.sh")"
boost_src_dir="$project_dir"'/src/3rd_party/boost_1_65_1'
boost_install_dir="$project_dir"'/src/3rd_party/boost'
if [ ! -d "$boost_src_dir" ]
then
    echo "boost src directory doesn't exist: $boost_src_dir"
    exit 1
fi
if [ -d "$boost_install_dir" ]
then
    echo "removing existing boost install directory: $boost_install_dir"
    rm -rf "$boost_install_dir"
fi

if [[ "$(hostname)" = "c2001" ]]
then
    echo "cluster build environment..."
    export CC="/usr/local/gcc-7.2.0/bin/gcc"
    export CXX="/usr/local/gcc-7.2.0/bin/g++"
    export MPI_C="/usr/local/mpich-3.1.3/bin/mpicc"
    export MPI_CXX="/usr/local/mpich-3.1.3/bin/mpicxx"
    export MPIEXEC="/usr/local/mpich-3.1.3/bin/mpiexec"
    export CONDA_PYTHON="/ifs/loni/faculty/dong/mcp/utils/anaconda2/bin/python"
    gcc_ver="$("$CC" -dumpversion)"
else
    gcc_ver="$(gcc -dumpversion)"
fi
if [[ "$gcc_ver" < "5.4.0" ]]
then
    echo "needs gcc version at least 5.4.0"
    exit 1
fi

cd "$boost_src_dir"
if [ ! -e "bootstrap.sh" ]
then
    echo "bootstrap.sh not found in boost src directory: $boost_src_dir"
    exit 1
fi

user_config="./tools/build/src/user-config.jam"
# this config does not build mpi or parallel graph, need using mpi in user_config
# c11 features 
echo "using gcc : $gcc_ver : $CXX : <cxxflags>-std=c++11 ;" > "$user_config"
# python build on cluster needs some help
if [[ "$(hostname)" = "c2001" ]]
then
    echo "using python : 2.7 : $CONDA_PYTHON ;" >> "$user_config"
fi
# configure
./bootstrap.sh --prefix="$boost_install_dir"
if [ ! -e "b2" ]
then
    echo "b2 not found in boost src directory: $boost_src_dir"
    exit 1
fi
# build. output log to install dir
./b2 install > "$boost_src_dir/b2_install.log"
mv "$boost_src_dir/b2_install.log" "$boost_install_dir/b2_install.log"
echo "build log written to $boost_install_dir/b2_install.log"
