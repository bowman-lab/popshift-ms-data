#!/bin/bash

pocket_sel=loos-pocket-sel-ca.txt
resect=1.0
model_tag=resect-$resect/tica-500-msm-2000-k-75-nframes-20
extracts=extracted_scores/12xsmina

for rep in {1..3}; do
    while read mol pdb_id; do
        if [[ $mol == apo ]]; then
            echo ignoring $mol
        else
            rstruct=xtals/$mol/$pdb_id.pdb
            echo $rstruct
            if [[ -f $rstruct ]]; then
                lig_res_name=$(python extract_oneval_json.py xtals/ligand-name-resname.json $mol)

                echo $mol : $rstruct : $lig_res_name
                if [[ $mol == 1-methylpyrrole ]]; then
                    python compare_confs_ref.py \
                        --ligand-sel "resid == 303" \
                        $rstruct \
                        $pocket_sel \
                        t4l-$rep/binding/$model_tag/receptor \
                        t4l-$rep/binding/$model_tag/$extracts/$mol.pickle 
                else
                    python compare_confs_ref.py \
                        --ligand-sel "resname == \"$lig_res_name\"" \
                        $rstruct \
                        $pocket_sel \
                        t4l-$rep/binding/$model_tag/receptor \
                        t4l-$rep/binding/$model_tag/$extracts/$mol.pickle 
                fi
            else
                echo no structs for $mol
            fi
        fi
    done < xtals/chosen-structs.txt
    
done
