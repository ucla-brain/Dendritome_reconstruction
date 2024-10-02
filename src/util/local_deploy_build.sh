#!/usr/bin/env bash

script_dir="$(dirname "${BASH_SOURCE[0]}")"
project_dir="$(bash "$script_dir/get_project_dir.sh")"
cd "$project_dir"
current_branch=$(git rev-parse --abbrev-ref HEAD)
echo "project directory: $project_dir. current branch: $current_branch"
build_dir="$project_dir/local_deploy_build"
if [[ "$1" = "clean-build" ]]
then
    echo "removing $build_dir"
    rm -rf "$build_dir"
fi
echo "build directory: $build_dir"

mkdir -p "$build_dir"
cd "$build_dir"
cmake -DCMAKE_BUILD_TYPE=Release -DMPI=FALSE -DPYTHON=TRUE -DPYTHON_VERSION=3 "$project_dir"'/src'
make
make install
exit 0