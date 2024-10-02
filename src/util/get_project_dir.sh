#!/usr/bin/env bash

# this script just echos project directory
project_dir="$(dirname $(dirname $(dirname $(readlink -f $0))))"
echo "$project_dir"