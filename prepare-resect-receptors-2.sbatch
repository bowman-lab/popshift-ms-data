#!/bin/bash
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p all
#SBATCH -o prepare-resect-receptors-2.out

find find t4l-2/ -name '*.pdb' -type f -path '*receptor/*' > resect-receptor-pdbs-2.txt
parallel -n 128 prepare_receptor -r {} -o {}qt :::: resect-receptor-pdbs-2.txt

