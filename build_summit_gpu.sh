#!/bin/bash

# Setup Summit modules
source summit_modules.sh

# Set paths and compiler flags
export RRTMGP_ROOT=${HOME}/codes/rte-rrtmgp/branches/improve-cloud-optics-example-gpu
export RTE_KERNELS=openacc
export GPTL_ROOT=${PROJWORK}/cli115/gptl
export NCHOME=$OLCF_NETCDF_ROOT
export NFHOME=$OLCF_NETCDF_FORTRAN_ROOT
export FC=mpif90
export FCFLAGS="-O2 -g -ta=nvidia,cc70,ptxinfo -Minfo=accel"

# Build RRTMGP library
cd $RRTMGP_ROOT/build
rm -f Makefile.conf
make clean
make -j8

# Build the rte-rrtmgp-clouds example
cd $RRTMGP_ROOT/examples/rte-rrtmgp-clouds
make clean
make -j8
