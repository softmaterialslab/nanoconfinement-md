#!/bin/bash
#SBATCH --nodes=1	
#SBATCH --ntasks-per-node=16
#SBATCH --time=0:30:00
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=IONS_LAMMPS
#SBATCH --output=out.log
#SBATCH --error=err.log

######  Module commands #####
module swap PrgEnv-intel PrgEnv-gnu
module load boost/gnu
module load gsl
module load lammps/gnu/7Aug19

######  Job commands go below this line #####
cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=16

#chmod 777 NAME

# -d refers to number of cores. this should match ntasks-per-node
# default is n = 4 and d = 16

time srun -n 1 -d 16 lmp_mpi < in.lammps
