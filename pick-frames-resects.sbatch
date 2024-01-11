#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o pick-frames-resects.out
#SBATCH -p all

frame_picker=~/raid/software/MSM-docking/pick_align_frames.py
for rep in {1..3}; do
    for k in 75; do
        for lag in 500; do
            for msm_lag in 1000 2000; do
                for resect in 0.2 0.4 0.6 0.8 1.0; do
                    inputdir=t4l-$rep
                    model=prot_masses.pdb
                    # change this if wanting to switch from using centers
                    dtraj=$inputdir/clustering/resect-$resect/t4l-$rep-backbone-all-chis-dihedrals-pocket-tica-lag-$lag-tica-reduced-k-$k-dtrajs.h5
                    msm=$inputdir/models/bb-all-chis-pocket-resect-$resect/tica-$lag-msm-$msm_lag/$k-clusters/mle/eq-probs.npy
                    mapping=$(dirname $msm)/mapping.json
                    if [[ -f $msm ]]; then
                        for nframes in 1 3 5 10 20; do
                            system_name=t4l-$rep/binding/resect-$resect/tica-$lag-msm-$msm_lag-k-$k-nframes-$nframes  # will also be outdir
                            traj_files=traj-list-$rep.txt
                            align_resids=pocket-resids-5.txt
                            selector_method=random
                            echo python $frame_picker --find-bonds 1.95,1.35\
                            --make-receptor-sel-chain-A\
                            --assignments $dtraj\
                            --align-resid-list $align_resids\
                            --number-frames $nframes\
                            --mapping $mapping\
                            $system_name\
                            $model\
                            $msm\
                            $selector_method\
                            '!hydrogen'\
                            $traj_files
                        done
                    fi
                done
            done
        done
    done
done > pick_frames_resects.txt
parallel -j 16 :::: pick_frames_resects.txt
