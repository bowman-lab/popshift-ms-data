from pathlib import Path
import pickle
import json
import numpy as np
import re


def extract_score_from_vina_pdbqts(pdbqts, search_re=re.compile(r'^REMARK VINA RESULT')):
    outlist = []
    for pdbqt in pdbqts:
        with pdbqt.open() as f:
            for line in f:
                if search_re.match(line):
                    outlist.append(float(line.split(maxsplit=4)[3]))
                    break
    return np.array(outlist, dtype='f4')


number_hi_scores = 10
top_p = Path('/project/bowmanlab/bnovak/ligand-binding/t4-lysozyme/dock/')
outdir_p = Path('for_meghan')
replica_rx = re.compile('replica.')
state_rx = re.compile('state[0-9]*')
regenerate_file_index = False
file_index_p = Path('dockdict-ligand-names.pickle')
target_runs = list(map(Path, (
   f'replica{i}_4000centers_python_scripting{suff}' for i in range(1,4)
    for suff in ('', '_additional_ligands', '_nbutyl')
)))

msm_p = Path('/project/bowmanlab/bnovak/ligand-binding/t4-lysozyme/MSM/')

ligns = set(j.name for i in target_runs for j in (top_p / i).iterdir())
dockdict = {}
if regenerate_file_index:
    for name in ligns:
        dockdict[name] = []
        print(name)
        for model_p in target_runs:
            print(model_p)
            ligdir = top_p / model_p / name
            if ligdir.is_dir():
                dockdict[name] += list(ligdir.rglob('replica0/*.pdbqt'))
            else:
                print(ligdir, 'not found, continuing.')


    with file_index_p.open('wb') as f:
        pickle.dump(dockdict, f)
else:
    with file_index_p.open('rb') as f:
        dockdict = pickle.load(f)

hi_scores = dict()

for name, ligpaths in dockdict.items():
    report_dir = outdir_p/name
    report_dir.mkdir(exist_ok=True)
    scores = extract_score_from_vina_pdbqts(ligpaths)
    sortinds = np.argsort(scores)
    winner_inds = sortinds[:number_hi_scores]
    hi_scores[name] = []
    for ix in winner_inds:
        ligpath = ligpaths[ix]
        ligpathstr = str(ligpath)
        rep = replica_rx.search(ligpathstr)[0]
        state = state_rx.search(ligpath.name)[0]
        # it will be unique, we only need first element of resultant generator
        statep = next((msm_p / (rep + '_1000')).rglob(f'centers/{state}.pdb'))
        outstatep = (report_dir/(rep+'-'+state)).with_suffix('.pdb')
        outlig = report_dir/(rep+'-'+ligpath.name)
        # copy pdbqt and state pdb to report directory
        outlig.write_text(ligpath.read_text())
        outstatep.write_text(statep.read_text())
        # save dictionary for results
        hi_scores[name].append(
            {
                'score': float(scores[ix]),
                'ligand_conf': str(outlig.relative_to(outdir_p)),
                'receptor_conf': str(outstatep.relative_to(outdir_p))
            }
        )
# save the high scores as the manifest to the conformations.
with (outdir_p/'manifest.json').open('w') as f:
    json.dump(hi_scores, f, indent=4)
