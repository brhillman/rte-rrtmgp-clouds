#!/bin/bash
#BSUB -P cli115
#BSUB -J rte-rrtmgp
#BSUB -eo rte-rrtmgp.%J
#BSUB -W 30
#BSUB -nnodes 1
#BSUB -I

# Parse arguments
# The "RUN" directory must be somewhere in /gpfs, or this will fail. jsrun cannot see anything but gpfs
project=cli115
RRTMGP_ROOT=${PWD}/../../ #${HOME}/codes/rte-rrtmgp/branches/improve-cloud-optics-example-gpu
RRTMGP_RUN=${PWD}/ref #/gpfs/alpine/scratch/${USER}/${project}/rrtmgp-clouds-gpu
for i in "$@"; do
    case $i in
        --rundir=*)
            RRTMGP_RUN="${i#*=}"
            shift
            ;;
        --rrtmgp-root=*)
            RRTMGP_ROOT="${i#*=}"
            shift
            ;;
        --project=*)
            project="${i#*=}"
            shift
            ;;
        *)
            echo "Argument $i not recognized."
            exit 1
            ;;
    esac
done


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
./rrtmgp_clouds \
rrtmgp-clouds.nc \
$RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc \
$RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc ${ncol} 

./rrtmgp_clouds \
rrtmgp-clouds.nc \
$RRTMGP_ROOT/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc \
$RRTMGP_ROOT/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc ${ncol}

echo "Done running in ${RRTMGP_RUN}"
