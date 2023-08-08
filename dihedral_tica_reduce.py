import os
from pathlib import Path
import pickle
import numpy as np
import jug


@jug.TaskGenerator
def mini_cluster(tica_filename, stride, n_clusters=50, overwrite=False):
    from enspara import ra
    from deeptime.clustering import KMeans

    cluster_filename = f"{os.path.splitext(tica_filename)[0]}-k-{n_clusters}-mini-cluster-object.pkl"
    if os.path.exists(cluster_filename) and not overwrite:
        print('cluster file already exists')
        return cluster_filename

    ragged_tica = ra.load(tica_filename)
    lengths = ragged_tica.lengths
    # note we need the astype call here because RAs always load sub-arrays as dtype objects.
    flatdata = np.concatenate(ragged_tica).astype(np.float32)
    # Try to save a little ram...
    del ragged_tica
    clusterer = KMeans(n_clusters)  # for future ref, additional params go here.
    # note this will correspond to the flattened tica data, and will therefore need to be reshaped by lengths.
    dtrajs = clusterer.fit_transform(flatdata)
    cm = clusterer.fetch_model()
    # since overwrite is broken, delete file and create new
    if os.path.exists(cluster_filename):
        os.remove(cluster_filename)
    with open(cluster_filename, 'wb') as f:
        pickle.dump(cm, f)

    cluster_dtraj_filename = f"{os.path.splitext(tica_filename)[0]}-k-{n_clusters}-dtrajs.h5"
    r = ra.RaggedArray(dtrajs, lengths=lengths)
    # since overwrite is broken, delete file and create new
    if os.path.exists(cluster_dtraj_filename):
        os.remove(cluster_dtraj_filename)
    ra.save(cluster_dtraj_filename, r)

    return cluster_filename


@jug.TaskGenerator
def tica_reduce(feature_filename, lag_time, tica_filename, data=None, var_cutoff=0.9,
                save_tica_obj=False, overwrite=False, resect_to_frac=None):
    from deeptime.decomposition import TICA
    from deeptime.util.types import to_dataset
    from enspara import ra
    # if output and we are not overriding, return filename
    if os.path.exists(tica_filename) and not overwrite:
        print('tica file already exists')
        return tica_filename
    
    # if data was passed in as an argument, don't read it from disk
    if data:
        pass
    elif resect_to_frac:
        data = []
        for i in feature_filename:
            feat_traj = np.load(i)
            length = feat_traj.shape[0]
            resect_length = int(length * resect_to_frac)
            data.append(feat_traj[:resect_length, :])
    else:
        data = [np.load(i) for i in feature_filename]

    # print(f'Data shape: {data[0].shape}')
    print(f'Number of trajectories: {len(data)}')
    # I assume lag is in frames
    fitter = TICA(lagtime=lag_time, var_cutoff=var_cutoff, scaling='commute_map')
    estimator = fitter.fit(data)
    tica_model = estimator.fetch_model()
    transformed_trajs = [tica_model.transform(trj) for trj in data]
    lengths = [len(trj) for trj in transformed_trajs]
    r = ra.RaggedArray(np.concatenate(transformed_trajs), lengths=lengths)
    del transformed_trajs
    ra.save(tica_filename, r)

    np.save(tica_filename.replace('tica-reduced.h5', 'tica-cumvar.npy'), tica_model.cumulative_kinetic_variance)

    # save out eigenvectors to get a sense of which features are being selected
    np.save(tica_filename.replace('tica-reduced.h5', 'tica-lsvs.npy'), tica_model.singular_vectors_left)
    np.save(tica_filename.replace('tica-reduced.h5', 'tica-rsvs.npy'), tica_model.singular_vectors_right)
    np.save(tica_filename.replace('tica-reduced.h5', 'tica-feat-corr.npy'), tica_model.feature_component_correlation)

    print('Number of dimensions saved is: ', tica_model.output_dimension, 'out of', data[0].shape[1])

    if save_tica_obj:
        with open(tica_filename.replace('tica-reduced.h5', 'tica-object.pkl'), 'wb') as f:
            pickle.dump(tica_model, f)
    return tica_filename


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

features_paths = {
    't4l-1': "recluster/t4l-1-backbone-all-dihedrals-pocket-feature-fns.txt",
    't4l-2': "recluster/t4l-2-backbone-all-dihedrals-pocket-feature-fns.txt",
    't4l-3': "recluster/t4l-3-backbone-all-dihedrals-pocket-feature-fns.txt"
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
        'features_list': features_paths[protein],
        'tica-lags':  [100, 200, 500, 1000],
        'recluster': False,
        'resects': list(np.arange(0.2, 1.2, 0.2).round(1))
        # 'k': [25, 50, 75, 100],
    }

outpre = '.'
var_cutoff = 0.9
feat_traj_suff = '.npy'
for protein, specs in to_cluster.items():
    features_list = Path(specs['features_list']).read_text().split()
    feat_trajs = [np.load(i) for i in features_list]
    for resect_to_frac in specs['resects']:
        overwrite = specs['recluster']
        OUT_STEM = outpre + f'/{protein}/clustering/resect-{resect_to_frac:.1f}'
        traj_paths = specs['traj_paths']
        top_path = specs['top_path']
        selstr = specs['selstr']
        stride = specs['stride']
        description = specs['description']
        data = []
        for feat_traj in feat_trajs:
            length = feat_traj.shape[0]
            resect_length = int(length * resect_to_frac)
            data.append(feat_traj[:resect_length, :]) 
   
        outstem_p = Path(OUT_STEM)
        if not outstem_p.is_dir():
            outstem_p.mkdir(parents=True)
            
        for lag_time in specs['tica-lags']:
            if 'chi_selstr' in specs.keys():
                if 'which_chis' in specs.keys():
                    tica_filename = f"{OUT_STEM}/{protein}-backbone-{specs['which_chis']}-chis-dihedrals-{description}-tica-lag-{lag_time}-tica-reduced.h5"
                else:
                    tica_filename = f"{OUT_STEM}/{protein}-backbone-chi1-dihedrals-{description}-lag-{lag_time}-tica-reduced.h5"
            else:
                tica_filename = f"{OUT_STEM}/{protein}-backbone-dihedrals-{description}-lag-{lag_time}-tica-reduced.h5"

            tica_filename = tica_reduce(features_list, lag_time, tica_filename, var_cutoff=var_cutoff, data=data, overwrite=overwrite, save_tica_obj=True, resect_to_frac=resect_to_frac)

            if 'k' in specs.keys():
                for k in specs['k']:
                    mini_cluster(tica_filename, stride, n_clusters=k, overwrite=overwrite)
