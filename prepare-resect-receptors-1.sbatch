#!/bin/bash
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p all
#SBATCH -o prepare-resect-receptors.out

find find t4l-1/ -name '*.pdb' -type f -path '*receptor/*' > resect-receptor-pdbs-1.txt
parallel -n 128 prepare_receptor -r {} -o {}qt :::: resect-receptor-pdbs-1.txt

