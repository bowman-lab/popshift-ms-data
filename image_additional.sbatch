#!/bin/bash
#SBATCH -p all
#SBATCH -n 1
#SBATCH -o reimage-trjconv.log

top=$(pwd)
for i in {1..3}; do
    for j in {1..5}; do
        cd replica$i/Production_MD$j
        echo "1 1" | gmx trjconv -f md2.*.xtc -o md_centered.xtc -n ../prev/merged.ndx -s md2.tpr -pbc mol -ur compact -center
        cd $top
    done
done


