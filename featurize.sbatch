#!/bin/bash
#SBATCH -p all
#SBATCH -N 8
#SBATCH --ntasks-per-node 8
#SBATCH -J t4l-df
#SBATCH -o log/dihe-feat-%A-%J.log

echo "executing task "
srun jug execute dihedral_tica_featurize.py
echo 'done with task '
