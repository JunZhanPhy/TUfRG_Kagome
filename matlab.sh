#!/bin/sh
#SBATCH -J 0.0_0.2
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -p regular

#module load mpi/intelmpi/2019u3
#module load mathlib/lapack/intel/3.4.2
ulimit -s unlimited
echo "Job starts at" date
matlab  -nodesktop -nosplash -nodisplay -r main
#~/soft/matlab_R2020a/bin/matlab  -nodesktop -nosplash -nodisplay -r main
#$EXE <solver.m>cal.out
echo "Job ends at" date
