#!/bin/bash

set -eu

main_dir=$(cd $(dirname $0); pwd -P)
src_dir="$main_dir/src"
build_dir="$main_dir/.build"
bin_dir="$main_dir/bin"

mkdir -p "$build_dir"
cd "$build_dir"
cmake ../
cmake --build .
cmake --install .

echo "---------------------------------------"
echo "Build and Install completed successfully!"
echo "Executable is in: $bin_dir"
echo "---------------------------------------"
