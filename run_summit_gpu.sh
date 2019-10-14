#!/bin/bash
#BSUB -P cli115
#BSUB -J rte-rrtmgp
#BSUB -eo rte-rrtmgp.%J
#BSUB -W 30
#BSUB -nnodes 1
#BSUB -I

# The "RUN" directory must be somewhere in /gpfs, or this will fail. jsrun cannot see anything but gpfs
project=cli115
RRTMGP_ROOT=${HOME}/codes/rte-rrtmgp/branches/improve-cloud-optics-example-gpu
RRTMGP_RUN=/gpfs/alpine/scratch/${USER}/${project}/rrtmgp-clouds-gpu

# Setup Summit modules
source summit_modules.sh

# Setup run directory in gpfs
mkdir -p $RRTMGP_RUN
cd $RRTMGP_RUN
rm *.nc
cp $RRTMGP_ROOT/examples/rte-rrtmgp-clouds/rrtmgp_clouds .
cp $RRTMGP_ROOT/examples/rte-rrtmgp-clouds/rrtmgp-clouds.nc .
cp $RRTMGP_ROOT/rrtmgp/data/*.nc .

# Run example problem with arbitrary number of copied columns
ncol=3500
echo "Run test for ncol = ${ncol}..."
jsrun -n 1 -c 1 -a 1 -g 1 nvprof --openmp-profiling off -o rrtmgp-clouds-lw.nvprof.%h.%p ./rrtmgp_clouds \
    rrtmgp-clouds.nc \
    $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc \
    $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc ${ncol} 
jsrun -n 1 -c 1 -a 1 -g 1 nvprof --openmp-profiling off -o rrtmgp-clouds-sw.nvprof.%h.%p ./rrtmgp_clouds \
    rrtmgp-clouds.nc \
    $RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc \
    $RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc ${ncol}

echo "Done running in ${RRTMGP_RUN}"
