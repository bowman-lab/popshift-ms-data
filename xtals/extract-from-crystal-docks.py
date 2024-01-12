import json
import numpy as np
from jug import TaskGenerator
from vina import Vina
import subprocess as sp
import re

from pathlib import Path


def extract_score_from_vina_pdbqt(pdbqts, search_re=re.compile(r'^REMARK VINA RESULT')):
    outlist = []
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[3]))
                    break
    return np.array(outlist, dtype='f4')

def extract_score_from_smina_pdbqt(pdbqts, search_re=re.compile(r'^REMARK minimizedAffinity')):
    outlist = []
    # print(pdbqts)
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[2]))
                    break
    return np.array(outlist, dtype='f4')


number_pattern = re.compile('\d+')

receptor_p = Path('n-butylbenzene/186L.pdbqt')
ligands_dir = Path('../all_ligands')
out_dir_p = Path('186L-docks')
box_center = (0, 0, 0)
box_size = (12, 12, 12)
exhaust = 32

out_dir_p.mkdir(exist_ok=True, parents=True)
outdict = {}
for ligp in ligands_dir.rglob('*.pdbqt'):
    output_pose_p = out_dir_p/ligp.name
    extracted_score = extract_score_from_smina_pdbqt([output_pose_p])
    outdict[ligp.stem] = float(extracted_score[0])

with (out_dir_p/'extracts.json').open('w') as f:
    json.dump(outdict, f)
