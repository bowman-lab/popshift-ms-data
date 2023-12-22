import numpy as np
from pathlib import Path
import pyemma
import os
import jug
import time

@jug.TaskGenerator
def cluster_subset(tica_filename, trj_len_cutoff, n_clusters, trial, output_stem, test_size=0.5, max_iter=500):
    from sklearn.model_selection import train_test_split
    from enspara import ra
    from deeptime.clustering import KMeans

    tica_trjs = ra.load(tica_filename)
    long_tica_trjs = []
    full_lengths = []
    for t in tica_trjs:
        length = len(t)
        if len(t) >= trj_len_cutoff:
            long_tica_trjs.append(t.astype(np.float32))
            full_lengths.append(length)

    del tica_trjs
    # Train test split
    train_tica_trjs, test_tica_trjs, train_lengths, test_lengths = train_test_split(long_tica_trjs, full_lengths, test_size=test_size)

    del long_tica_trjs

    # time clustering
    start = time.time()

    clusterer = KMeans(n_clusters=n_clusters, max_iter=max_iter)
    # note because we have to concat this will come off as one uniform hunk, not an ra
    flat_train_assigns = clusterer.fit_transform(np.concatenate(train_tica_trjs))
    end = time.time()
    print(f'clustering took {end - start} seconds with k={n_clusters}.')
    train_assignments = ra.RaggedArray(flat_train_assigns, lengths=train_lengths)
    train_assignment_filename = f'{output_stem}-k-{n_clusters}-split-{trial}-train-ra.h5'

    # make directory if it does not already exist
    os.makedirs(os.path.dirname(train_assignment_filename), exist_ok=True)

    # np.save(train_assignment_filename, train_assignments)
    ra.save(train_assignment_filename, train_assignments)

    flat_test_assigns = clusterer.transform(np.concatenate(test_tica_trjs))
    test_assignment_filename = f'{output_stem}-k-{n_clusters}-split-{trial}-test-ra.h5'
    # np.save(test_assignment_filename, test_assignments)
    ra.save(test_assignment_filename, ra.RaggedArray(flat_test_assigns, lengths=test_lengths))

    return train_assignment_filename, test_assignment_filename

@jug.TaskGenerator
def vamp2_score(assignment_filenames, lag_time):
    from enspara import ra
    train_assignment_filename, test_assignment_filename = assignment_filenames
    dtrajs_train = list(ra.load(train_assignment_filename))
    dtrajs_test = list(ra.load(test_assignment_filename))

    # VAMP-2
    # time clustering
    try:
        start = time.time()
        pyemma_msm = pyemma.msm.estimate_markov_model(dtrajs_train, lag=lag_time, score_method='VAMP2', score_k=10)
        end = time.time()
        print(f'MSM fitting took {end-start} for {os.path.basename(train_assignment_filename)}')

        start = time.time()
        vamp2_train_score = pyemma_msm.score(dtrajs_train, score_method='VAMP2', score_k=10)
        end = time.time()
        print(f'VAMP scoring took {end-start} for {os.path.basename(train_assignment_filename)}')

        vamp2_test_score = pyemma_msm.score(dtrajs_test, score_method='VAMP2', score_k=10)

        return vamp2_train_score, vamp2_test_score
    except:
        print(f'trial failed for {train_assignment_filename}')
        return np.nan, np.nan

@jug.TaskGenerator
def save_vamp2_scores(vamp2_scores, output_basename):
    print(vamp2_scores)
    vamp2_train_scores = [[s[0] for s in scores_for_k] for scores_for_k in vamp2_scores]
    vamp2_test_scores = [[s[1] for s in scores_for_k] for scores_for_k in vamp2_scores]
    np.save(f'{output_basename}-train.npy', vamp2_train_scores)
    np.save(f'{output_basename}-test.npy', vamp2_test_scores)
    return None



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
        'output_dir': f'{protein}/clustering/vamp-scan',
        'clustering_pre': f'{protein}/clustering/resect-1.0',
        'description': 'pocket',
        'features_list': features_paths[protein],
        'tica-lags':  [100, 200, 500, 1000],
        'k': [25, 50, 75, 100, 200, 300, 400, 500],
    }

n_trials = 10
trj_len_cutoff = 10 * 50  # in frames (50 ns * 50 frames / ns)


for protein, specs in to_cluster.items():
    clustering_pre = specs['clustering_pre']
    output_dir = specs['output_dir']
    outstem_p = Path(output_dir)
    description = specs['description']
    if not outstem_p.is_dir():
        outstem_p.mkdir(parents=True)
    for lag_time in specs['tica-lags']:
        if 'chi_selstr' in specs.keys():
            if 'which_chis' in specs.keys():
                tica_filename = f"{clustering_pre}/{protein}-backbone-{specs['which_chis']}-chis-dihedrals-{description}-tica-lag-{lag_time}-tica-reduced.h5"
            else:
                tica_filename = f"{clustering_pre}/{protein}-backbone-chi1-dihedrals-{description}-lag-{lag_time}-tica-reduced.h5"
        else:
            tica_filename = f"{clustering_pre}/{protein}-backbone-dihedrals-{description}-lag-{lag_time}-tica-reduced.h5"
        print(Path(tica_filename).is_file(), tica_filename)
        output_stem = f"{output_dir}/{os.path.basename(tica_filename).split('.')[0]}"
        print(output_stem)
        vamp2_scores = []

        for k in specs['k']:
            scores_for_k = []
            for trial in range(n_trials):
                assignment_filenames = cluster_subset(
                    tica_filename, trj_len_cutoff, k, trial, output_stem)
                scores_for_k.append(vamp2_score(assignment_filenames, lag_time))
            vamp2_scores.append(scores_for_k)

        scan_description = '-'.join(str(k) for k in specs['k'])
        output_basename = f"{output_stem}-vamp-scan-lagtime-{lag_time}-k-{scan_description}"
        save_vamp2_scores(vamp2_scores, output_basename)

