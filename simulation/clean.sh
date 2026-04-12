#!/bin/sh

set -e

main_dir=$(cd $(dirname $0); pwd -P)
build_dir=$main_dir/.build
bin_dir=$main_dir/bin

rm -rf $build_dir $bin_dir
