#!/bin/bash

script_dir="$(dirname "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "../../util/get_project_dir.sh")"
cd "$project_dir/benchmark/cluster_runtime"
qsub -cwd "cv_vs_eigen.sh"
