#!/bin/bash

script_dir="$(dirname "BASH_SOURCE[0]")"
cd "$script_dir"
project_dir="$(bash "../util/get_project_dir.sh")"
cd "$project_dir/grid_build_master/parallel"
qsub -cwd "par_cc_3d.sh"