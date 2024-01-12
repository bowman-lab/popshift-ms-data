#!/bin/bash

lines=1000
# convert_one(){
#     pdbqt_basename=$1
#     if [[ ! -f $pdbqt_basename.pdb ]]; then
#         obabel -ipdbqt -opdb $pdbqt_basename.pdbqt -O$pdbqt_basename.pdb
#     fi
# }
convert_several(){
    for pdbqt_basename in "$@"; do
        obabel -ipdbqt -opdb $pdbqt_basename.pdbqt -O$pdbqt_basename.pdb
    done
}
# export -f convert_one
export -f convert_several
runname=12xsmina
find t4l-?/binding/resect-*/*/$runname/ -name '*.pdbqt' > potential_converts.txt
parallel -l $lines --bar -j 16 convert_several {.} :::: potential_converts.txt