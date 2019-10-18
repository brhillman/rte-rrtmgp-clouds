#!/bin/bash

# Setup Summit modules
source summit_modules.sh

# Parse arguments
gpu=0
for i in "$@"; do
    case $i in
        --gpu=*)
            gpu="${i#*=}"
            shift
            ;;
        --help)
            echo "usage: `basename $0` <args>"
            echo ""
            echo "Arguments:"
            echo ""
            echo "--gpu=<1, 0>          Flag for whether or not to build for GPU (0: CPU, 1: GPU)"
            echo "--help                Display help message"
            echo ""
            echo "Author: Ben Hillman (bhillma@sandia.gov)"
            exit 0
            ;;
        *)
            echo "Argument $i not recognized."
            exit 1
            ;;
    esac
done



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
