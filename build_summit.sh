#!/bin/bash

# Setup Summit modules
source summit_modules.sh

gpu=0
if [ $# -eq 1 ]; then
    gpu=$1
elif [ $# -gt 1]; then
    echo "usage: `basename $0` <GPU>"
    echo "       GPU = [0: use CPU; 1: use GPU]"
    exit 1
fi

# Set paths and compiler flags
export RRTMGP_ROOT=${HOME}/codes/rte-rrtmgp/branches/improve-cloud-optics-example-gpu
export FC=mpif90
if [ ${gpu} -eq 1 ]; then
    export RTE_KERNELS=openacc
    export FCFLAGS="-O2 -g -ta=nvidia,cc70,ptxinfo -Minfo=accel"
else
    unset RTE_KERNELS
    export FCFLAGS="-O2 -g"
fi
export GPTL_ROOT=${PROJWORK}/cli115/gptl
export NCHOME=$OLCF_NETCDF_ROOT
export NFHOME=$OLCF_NETCDF_FORTRAN_ROOT

# Build RRTMGP library
cd $RRTMGP_ROOT/build
rm -f Makefile.conf
make clean
make -j8

# Build the rte-rrtmgp-clouds example
cd $RRTMGP_ROOT/examples/rte-rrtmgp-clouds
make clean
make -j8
