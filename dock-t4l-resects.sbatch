#!/bin/bash
#SBATCH --array=0-1000
#SBATCH -n 1
#SBATCH -p all
#SBATCH -J lysos
#SBATCH -o dockouts/resect-%A-%a.out
popshift_path=~/raid/software/PopShift
for rep in {1..3}; do
	inputdir=t4l-$rep
	for lag in 500; do
		for msm_lag in 2000; do
			for k in 75; do
				for resect in 1.0; do
					for mtype in mle; do
						msm=$inputdir/models/bb-all-chis-pocket-resect-$resect/tica-$lag-msm-$msm_lag/$k-clusters/$mtype/eq-probs.npy
						ls $msm
						if [[ -f $msm ]]; then
							for nframes in 10; do
								receptor_dir=t4l-$rep/binding/resect-$resect/tica-$lag-msm-$msm_lag-k-$k-nframes-$nframes
								ls $receptor_dir
								echo running $receptor_dir
								jug execute -- $popshift_path/docking_parallel.py\
								 -d smina\
								 -e 32\
								 $receptor_dir/receptor\
								 $receptor_dir/'12xsmina'\
								 '0,0,0'\
								 '12,12,12'\
								 expanded-ligands/*/*.pdbqt
								echo finished $receptor_dir
					
							done

						fi
					done
				done
			done
		done
	done
done

