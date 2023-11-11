#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH -A r00458
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=fanbsun@iu.edu
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --job-name=10_1_-1_0.2_0.2_0.1_-0.01

# below are the modules you will need to compile the code on bigred200 (see README)
######  Module commands for inhouse code and lammps #####
module load lammps/2Aug2023


cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

# the following is a test simulation to check things are working, it simulates 0.1 M of +1 and -1 ions of diameter 0.714 nm confined within 3 nm
# -d refers to number of cores. this should match ppn in Line 2.
time srun -n 96 -d 1 lmp_mpi < in.lammps