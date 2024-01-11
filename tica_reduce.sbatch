#!/bin/bash
#SBATCH -p all
#SBATCH -N 4
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32GB
#SBATCH -c 8
#SBATCH -J t4l-clus
#SBATCH -o log/tica_reduce-%J-%A.log
echo starting jug task
srun jug execute dihedral_tica_reduce.py
echo 'done with task '
