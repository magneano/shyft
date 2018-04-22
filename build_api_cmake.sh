#!/bin/bash

set -e
exec 3>&1 4>&2
SHYFT_WORKSPACE=${SHYFT_WORKSPACE:=$(readlink --canonicalize --no-newline `dirname ${0}`/..)}
SHYFT_DEPENDENCIES_DIR=${SHYFT_DEPENDENCIES_DIR:=${SHYFT_WORKSPACE}/shyft_dependencies}
bash build_support/build_dependencies.sh
export LD_LIBRARY_PATH=${SHYFT_DEPENDENCIES_DIR}/lib:$LD_LIBRARY_PATH
export PATH=${SHYFT_WORKSPACE}/miniconda/bin:$PATH
mkdir -p build
cd build
cmake ..
make -j 4
make install
