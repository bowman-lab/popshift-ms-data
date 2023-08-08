import os
from glob import glob

import jug
import mdtraj as md
import numpy as np
import pyemma
from pathlib import Path


@jug.TaskGenerator
def featurize_write(featurizer, traj_name, stride, traj_feature_suffix='.h5', parent=None, 
                    rename_depth=-4, chunksize=2000):
    features = pyemma.coordinates.load(traj_name, features=featurizer, stride=stride, chunksize=chunksize)
    if parent:
        tp = Path(traj_name).with_suffix(traj_feature_suffix)
        traj_specific_feature_p = (parent / Path(*tp.parts[rename_depth:])).resolve()
        traj_specific_feature_fn = str(traj_specific_feature_p)
        traj_specific_feature_p.parent.mkdir(exist_ok=True, parents=True)
    else:
        traj_specific_feature_fn = traj_name.replace('.xtc', traj_feature_suffix)
    print(traj_specific_feature_fn)
    np.save(traj_specific_feature_fn, features)
    return traj_specific_feature_fn



@jug.TaskGenerator
def write_feature_traj_fns(feature_traj_list_name, features_list):
    with open(feature_traj_list_name, 'w') as f:
        f.write('\n'.join(features_list))


topologies = {
    't4l-1': 'prot_masses.pdb',
    't4l-2': 'prot_masses.pdb',
    't4l-3': 'prot_masses.pdb',
}

trajectory_paths = {
    't4l-1': [
        'traj-list-1.txt'
    ],
    't4l-2': [
        'traj-list-2.txt'
    ],
    't4l-3': [
        'traj-list-3.txt'
    ]
}

pocket_resids = np.loadtxt('pocket-resids-5.txt', dtype=int)

to_cluster = {}
for protein in trajectory_paths.keys():
    print(protein)
    to_cluster[protein] = {
        'traj_paths': trajectory_paths[protein],
        'top_path': topologies[protein],
        'stride': 1,
        'selstr': ' or '.join(f'residue {r}' for r in pocket_resids),
        'chi_selstr': ' or '.join(f'residue {r}' for r in pocket_resids),
        'which_chis': "all",
        'description': 'pocket',
        'tica-lags': [100, 200, 500, 1000],
        'k': [25, 50, 75, 100, 150, 200, 300, 400, 500],
    }

var_cutoff = 0.9
feat_traj_suff = '.npy'
overwrite = True
OUT_STEM = "recluster"
for protein, specs in to_cluster.items():
    traj_paths = specs['traj_paths']
    top_path = specs['top_path']
    selstr = specs['selstr']
    stride = specs['stride']
    description = specs['description']

    if 'chi_selstr' in specs.keys():
        include_sidechains = True
        chi_selstr = specs['chi_selstr']
        output_filename = f"{OUT_STEM}/{protein}-backbone-{specs['which_chis']}-dihedrals-{description}-feature-fns.txt"
        if 'which_chis' in specs.keys():
            which_chis = specs['which_chis']
        else:
            which_chis = 'chi1'
    else:
        chi_selstr = selstr
        output_filename = f'{OUT_STEM}/{protein}-backbone-dihedrals-{description}-feature-fns.txt'
        include_sidechains = False

    # If feature filename already exists, continue to next component of pipeline.
    if os.path.exists(output_filename):
        if overwrite:
            print('Features filename list found, but overwrite set to "True", so overwriting.')
        else:
            print('Feature file already exists; refusing to redo.')
            exit(1) 
    pdb = md.load(top_path)
    feat = pyemma.coordinates.featurizer(pdb)
    feat.add_backbone_torsions(selstr=selstr, cossin=True, periodic=False)
    if include_sidechains:
        feat.add_sidechain_torsions(selstr=chi_selstr, cossin=True, periodic=False, which=which_chis)

    # save out description of features
    os.makedirs(OUT_STEM, exist_ok=True)
    np.save(output_filename.replace('feature-fns.txt', 'feature-descriptions.npy'), feat.describe())

    if Path(traj_paths[0]).suffix == '.txt':
        traj_list = []
        for traj_list_path in traj_paths:
            traj_list += Path(traj_list_path).read_text().strip().split()
    else:  # assume these are trajectory filenames
        traj_list = sorted(list(np.concatenate([glob(traj_path) for traj_path in traj_paths])))

    print('Beginning to generate tasks for reading, featurizing, and writing per trajectory.')
    # spawn a bunch of jug tasks to read and featurize individual trajectories, saving them to independent h5 files
    # NOTE: if parent is defined, will write to new directory structure rename_depth below top levels of traj files
    # within parent. If not, it will write a 'feature.h5' in the same path as the parent of the trajectory file.
    features_list = [featurize_write(feat, traj_name, stride, traj_feature_suffix=feat_traj_suff, parent=OUT_STEM)
                     for traj_name in traj_list]
    # write the file names for each feature traj to a list that can be read by subsequent (clustering) scripts.
    write_feature_traj_fns(output_filename, features_list) 
