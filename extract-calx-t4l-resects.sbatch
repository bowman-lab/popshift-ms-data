#!/bin/bash
#SBATCH -N 1
#SBATCH -c 8 
#SBATCH -n 8
#SBATCH -p all
#SBATCH -o extract-calx-t4l-resects-all.out

msmbindpath=~/raid/software/MSM-docking
calx=$msmbindpath/msm_binder.py
extract=$msmbindpath/extract-scores.py
ncores=8 # this'll be the number of cores given to each job

for rep in {1..3}; do
	top=t4l-$rep
	binding=$top/binding
	models=$top/models
	for k in 75; do
		for lag in 500; do
			for msm_lag in 1000 2000; do
				for resect in 0.2 0.4 0.6 0.8 1.0; do
					for nframes in 1 3 5 10 20; do
						dock_results=$binding/resect-$resect/tica-$lag-msm-$msm_lag-k-$k-nframes-$nframes
						for runname in 12xsmina; do
							extracted_scores=$dock_results/extracted_scores/$runname
							dock_run=$dock_results/$runname
							echo starting $dock_run
							docktype=smina
							python $extract -n $ncores -t $docktype $dock_run
							for modeltype in mle; do
								model_eq=$models/bb-all-chis-pocket-resect-$resect/tica-$lag-msm-$msm_lag/$k-clusters/$modeltype/eq-probs.npy
								if ls $extracted_scores/*.h5; then
									python $calx -n $ncores --reweighted-eq \
										--out $extracted_scores/binding-calx-$modeltype \
										--reweighted-eq 1 \
										bin-samples $model_eq \
										$extracted_scores/*.h5 &
								else
									echo $(readlink -f $dock_results/$runname) : extraction failed.
								fi
							done
						done
					done
				done
				wait
			done
		done
	done
done
