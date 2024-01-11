#!/bin/bash
#SBATCH -p all
#SBATCH -N 14
#SBATCH --ntasks-per-node 1
#SBATCH -c 48
#SBATCH -J t4l-vamp
#SBATCH -o log/vamp-%A.log
#SBATCH -e log/vamp-%A.log

echo "executing task "
srun jug execute parallel_vamp_scan_train_test_split.py
echo 'done with task '
