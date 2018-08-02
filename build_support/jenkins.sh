#!/bin/bash

# This script will attempt to:
#  build shyft,
#  run tests,
#  generate coverage reports
#  make a conda package
#  upload the package using anaconda
#
# The script may be invoked by build robots such as Jenkins.

# Set any host-specific configuration
host=$(hostname)
case $host in
    'D40060-ubuntu')
	export WORKSPACE=/home/jenkins/workspace;
	export PATH=/home/u40420/projects/cmake/bin:$PATH
	;;
    'oslxpsht002p.energycorp.com')
	export WORKSPACE=/var/lib/jenkins/workspace;
	;;
    *) echo "Error, no configuration defined for this host (${host})";
       exit 1
       ;;
esac;


# Make sure the script exits if any command fails
set -e

# Set paths and activate conda environment
source $WORKSPACE/miniconda/etc/profile.d/conda.sh
conda activate shyft_env
export SHYFT_DEPENDENCIES_DIR=$WORKSPACE/shyft_dependencies
export PYTHONPATH=$WORKSPACE/shyft

# Build shyft
cd $WORKSPACE/shyft
mkdir -p build
cd build
cmake ..
make -j 4 CMAKE_VERBOSE_MAKEFILE=0

# Run tests
make install
make test
cd ..
python -c "import shyft; shyft.print_versions()"
export SHYFT_SKIP_OPENDAP_TEST=1
nosetests shyft/tests --with-coverage --cover-package=shyft.repository --cover-package=shyft.orchestration --cover-package=shyft.api

# Convert coverage report to xml to be parsed by Cobertura plugin
coverage xml

# Build conda package
numpy_version=$(python -c "import numpy; print(numpy.version.short_version)")
shyft_minor=$(git rev-list --count HEAD)
export SHYFT_VERSION="4.6.${shyft_minor}"
conda build --numpy $numpy_version conda_recipe

# Upload conda package
filename=$(conda build --numpy $numpy_version --output conda_recipe)
anaconda upload --force --user statkraft ${filename}
