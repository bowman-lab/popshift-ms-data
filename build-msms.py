from pathlib import Path
from functools import partial
import numpy as np
from enspara import ra
from enspara.msm import builders, msm
import deeptime

chosen_k = [75]
lags = [500]
msm_lags = [2000]
reps = range(1,4)

for rep in reps:
    for lag in lags:
        for k in chosen_k:
            basep = Path(f't4l-{rep}')
            modelp = basep/'models'/'bb-all-chis-pocket'/f'tica-msm-{lag}'/f'{k}-clusters'
            if not modelp.is_dir():
                modelp.mkdir(parents=True, exist_ok=True)
            assignp = basep/'clustering'/f't4l-{rep}-backbone-all-chis-dihedrals-pocket-tica-lag-{lag}-tica-reduced-k-{k}-dtrajs.h5'
            if not assignp.is_file():
                print(assignp, ': not found. Skipping.')
                continue
            assigns = ra.load(str(assignp))
            mlefitter = deeptime.markov.msm.MaximumLikelihoodMSM(lagtime=lag, connectivity_threshold=1/k)
            mle = mlefitter.fit_fetch([a.astype(int) for a in assigns])
            print(mle.n_connected_msms)
            mle_tprobs = mle.transition_matrix
            mle_eq_probs = mle.stationary_distribution
            mledir = modelp/'mle'
            mledir.mkdir(exist_ok=True)
            np.save(str(mledir/'tprobs.npy'), mle_tprobs)
            np.save(str(mledir/'eq-probs.npy'), mle_eq_probs)
            normsm = msm.MSM(lag_time=lag, trim=True,
                             method=partial(builders.normalize, prior_counts=1/k, calculate_eq_probs=True))
            normsm.fit(assigns)
            normsmp = modelp/'norm'
            # if normsmp.is_dir():
            #     for fp in normsmp.iterdir():
            #         fp.unlink()
            #     normsmp.rmdir()
            normsm.save(str(normsmp), force=True)
            # get a version of eq probs that are numpy binary format at same overall path.
            eq_probs = np.loadtxt(str(normsmp/'eq-probs.dat'))
            np.save(str(normsmp/'eq-probs.npy'), eq_probs)
