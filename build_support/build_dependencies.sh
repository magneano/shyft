#!/bin/bash
export SHYFT_WORKSPACE=${SHYFT_WORKSPACE:=$(readlink --canonicalize --no-newline `dirname ${0}`/../..)}
# to align the cmake support:
SHYFT_DEPENDENCIES_DIR=${SHYFT_DEPENDENCIES_DIR:=${SHYFT_WORKSPACE}/shyft_dependencies}
armadillo_name=armadillo-9.200.6
dlib_name=dlib-19.16
boost_ver=1_68_0
pybind11_ver=v2.2.3
miniconda_ver=4.5.4
cmake_common="-DCMAKE_INSTALL_MESSAGE=NEVER"
echo ---------------
echo Update/build shyft dependencies
echo SHYFT_WORKSPACE........: ${SHYFT_WORKSPACE}
echo SHYFT_DEPENDENCIES_DIR.: ${SHYFT_DEPENDENCIES_DIR}
echo PACKAGES...............: miniconda ${miniconda_ver} w/shyft_env, doctest, boost_${boost_ver}, ${armadillo_name}, ${dlib_name} 

# A helper function to compare versions
function version { echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }'; }

# the current versions we are building
mkdir -p ${SHYFT_DEPENDENCIES_DIR}
cd ${SHYFT_DEPENDENCIES_DIR}

if [ ! -d ${armadillo_name} ]; then 
    echo Building ${armadillo_name}
    if [ ! -f ${armadillo_name}.tar.xz ]; then 
        wget  http://sourceforge.net/projects/arma/files/${armadillo_name}.tar.xz
    fi;
    tar -xf ${armadillo_name}.tar.xz
    pushd ${armadillo_name}
    cmake . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DDETECT_HDF5=false -DCMAKE_INSTALL_LIBDIR=lib ${cmake_common}
    make install
    popd
fi;
echo Done ${armadillo_name}

if [ ! -d ${dlib_name} ]; then
    echo Building ${dlib_name}
    if [ ! -f ${dlib_name}.tar.bz2 ]; then
        wget http://dlib.net/files/${dlib_name}.tar.bz2
    fi;
    tar -xf ${dlib_name}.tar.bz2
    pushd ${dlib_name}
    mkdir -p build
    dlib_cfg="-DDLIB_PNG_SUPPORT=0 -DDLIB_GIF_SUPPORT=0 -DDLIB_LINK_WITH_SQLITE3=0 -DDLIB_NO_GUI_SUPPORT=1 -DDLIB_DISABLE_ASSERTS=1 -DDLIB_JPEG_SUPPORT=0 -DDLIB_USE_BLAS=0 -DDLIB_USE_LAPACK=0 -DBUILD_SHARED_LIBS=ON"
    cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DCMAKE_INSTALL_LIBDIR=lib ${cmake_common} ${dlib_cfg} && cmake --build . --target install
    popd
fi;
echo Done ${dlib_name}

if [ ! -d doctest ]; then
    echo Building doctest
    git clone https://github.com/onqtam/doctest
    pushd doctest
    cmake . -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} ${cmake_common}
    make install
    popd
fi;
echo Done doctest

cd ${SHYFT_WORKSPACE}
# we cache minconda on travis, so the directory exists, have to check for bin dir:
if [ ! ${CONDA_PREFIX} ]; then
if [ ! -d miniconda/bin ]; then
    echo Building miniconda
    if [ -d miniconda ]; then
        rm -rf miniconda
    fi;
    if [ ! -f miniconda.sh ]; then
        wget  -O miniconda.sh http://repo.continuum.io/miniconda/Miniconda3-${miniconda_ver}-Linux-x86_64.sh
    fi;
    bash miniconda.sh -b -p ${SHYFT_WORKSPACE}/miniconda

    # Update conda to latest version, assume we start with 4.3 which
    # requires PATH to be set
    OLDPATH=${PATH}
    export PATH="${SHYFT_WORKSPACE}/miniconda/bin:$PATH"

    old_conda_version=$(conda --version | sed "s/conda \(.*\)/\1/")
    echo "Old conda version is ${old_conda_version}"
    if [[ $(version ${old_conda_version}) -ge $(version "4.4") ]]; then
	PATH=$OLDPATH
	source ${SHYFT_WORKSPACE}/miniconda/etc/profile.d/conda.sh
	conda activate
    else
	source activate
    fi
    conda config --set always_yes yes --set changeps1 no
    conda update conda

    new_conda_version=$(conda --version | sed "s/conda \(.*\)/\1/")
    echo "New conda version is ${new_conda_version}"
    if [[ $(version ${old_conda_version}) -lt $(version "4.4") &&
	      $(version ${new_conda_version}) -ge $(version "4.4") ]]; then
	PATH=$OLDPATH
	source ${SHYFT_WORKSPACE}/miniconda/etc/profile.d/conda.sh
	conda activate
    fi

    conda install numpy
    conda create -n shyft_env python=3.6 pyyaml numpy netcdf4 cftime gdal matplotlib requests nose coverage pip shapely  pyproj
    ln -s ${SHYFT_WORKSPACE}/miniconda/include/python3.6m ${SHYFT_WORKSPACE}/miniconda/include/python3.6
    ln -s ${SHYFT_WORKSPACE}/miniconda/envs/shyft_env/include/python3.6m ${SHYFT_WORKSPACE}/miniconda/envs/shyft_env/include/python3.6 
fi;
echo Done minconda
export PATH="${SHYFT_WORKSPACE}/miniconda/bin:$PATH"
fi;



cd ${SHYFT_DEPENDENCIES_DIR}
if [ ! -d boost_${boost_ver} ]; then
    echo Building boost_${boost_ver}
    if [ ! -f boost_${boost_ver}.tar.gz ]; then
        wget -O boost_${boost_ver}.tar.gz http://sourceforge.net/projects/boost/files/boost/${boost_ver//_/.}/boost_${boost_ver}.tar.gz
    fi;
    tar -xf boost_${boost_ver}.tar.gz
    pushd boost_${boost_ver}
    ./bootstrap.sh --prefix=${SHYFT_DEPENDENCIES_DIR}
    boost_packages="--with-system --with-filesystem --with-date_time --with-python --with-serialization"
    ./b2 -j2 -d0 link=shared variant=release threading=multi ${boost_packages}
    ./b2 -j2 -d0 install threading=multi link=shared ${boost_packages}
    popd
fi;
echo  Done boost_${boost_ver}

cd ${SHYFT_DEPENDENCIES_DIR}
if [ ! -d pybind11 ]; then
    git clone https://github.com/pybind/pybind11.git
    pushd pybind11
    git checkout master
    git pull
    git checkout ${pybind11_ver} > /dev/null
    mkdir -p build
    cd build && cmake .. -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DPYBIND11_TEST=0 ${cmake_common} .. && cmake -P cmake_install.cmake
    popd
fi;

echo Done pybind11

cd ${SHYFT_WORKSPACE}
if [ -d shyft-data ]; then 
    pushd shyft-data
    git pull >/dev/null
    popd
else 
    git clone https://github.com/statkraft/shyft-data
fi;
echo Done shyft-data
#echo Update shyft/shyft/lib with all 3rd party .so so that rpath will work for python extensions
#mkdir -p ${SHYFT_WORKSPACE}/shyft/shyft/lib
#install  --preserve-timestamps --target=${SHYFT_WORKSPACE}/shyft/shyft/lib ${SHYFT_DEPENDENCIES_DIR}/lib/*.so.*

