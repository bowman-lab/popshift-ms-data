#!/bin/bash
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p all
#SBATCH -o prepare-resect-receptors.out

for i in {1..3}; do 
	for resect in 0.2 0.4 0.6 0.8 1.0; do 
		parallel -j 128 "prepare_receptor -r {} -o {}qt" ::: t4l-$i/binding/resect-$resect/*/*/*/*.pdb
	done
done

