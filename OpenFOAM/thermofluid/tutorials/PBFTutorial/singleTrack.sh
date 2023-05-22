#!/bin/sh
# All the information about queues can be obtained using 'sinfo'

# Slurm flags
#SBATCH -N 4
#SBATCH --ntasks-per-node 20
#SBATCH -t 72:00:00
#SBATCH --job-name=1pt5mps_200W

# Mail me on job start & end
####SBATCH --mail-user=gowthaman.parivendhan@ucdconnect.ie
###SBATCH --mail-type=BEGIN,END

cd $SLURM_SUBMIT_DIR

# Load the module for the specific version
module load openmpi
source ~/foam/foam-extend-4.0/etc/bashrc

# Run simulation in parallel
mpirun -n 80 buoyantLaserMeltFoam -parallel

# Reconstruct the files
reconstructPar
