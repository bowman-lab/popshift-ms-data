#!/bin/bash
#SBATCH -p all
#SBATCH -N 8
#SBATCH --ntasks-per-node=1
#SBATCH -c 32
#SBATCH -J t4l-clus
#SBATCH -o log/cluster-%J.log

echo 'start task '
srun jug execute dihedral_tica_clustering.py
echo 'done with task '
