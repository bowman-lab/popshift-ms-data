import loos
from loos import pyloos
from loos import torsion
from pathlib import Path
import pickle
import numpy as np
from pathlib import Path

tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
model_type = 'mle'
runtag = '12xsmina'
binding_name = Path(f'binding/resect-{resect}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
extracts = binding_name / f'extracted_scores/{runtag}'
calx_name = extracts / f'binding-calx-{model_type}'
sel = 'resid == 111'
analysis_tag = 'val111_chi1'
dihe_atom_names = ["N", "CA", "CB", "CG1"]
# dihe_atom_names.reverse()
mol = 'benzene'
for rep in range(1,4):
    rep_dir = Path(f't4l-{rep}') 
    receptor_dir = rep_dir / binding_name / 'receptor'
    out_dir = receptor_dir.parent/'receptor-observables'
    if not out_dir.is_dir():
        out_dir.mkdir(parents=True)
    sample_confs = rep_dir /extracts/f'{mol}.pickle'
    with open(sample_confs, 'rb') as f:
        ligand_poseps = pickle.load(f)

    raw_confs = {x.stem: x for x in receptor_dir.rglob('*.pdb')}
    # flatten, but in same order as flattened RA.
    sample_conf_list = [raw_confs[lp.stem] for lpl in ligand_poseps for lp in lpl]
    model = loos.createSystem(str(sample_conf_list[0]))
    val111 = loos.selectAtoms(model, sel)
    dihe_ags = [loos.selectAtoms(val111, f'name == "{dn}"') for dn in dihe_atom_names]
    vtraj = pyloos.VirtualTrajectory(
        *(pyloos.Trajectory(str(fn), model, subset=sel) for fn in sample_conf_list)
    )
    dihes = []
    for _ in vtraj:
        dihes.append(torsion(*map(lambda x: x[0], dihe_ags)))
    np.save(str(out_dir/f'{analysis_tag}.npy'), dihes)
