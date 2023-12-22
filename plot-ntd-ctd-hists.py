import numpy as np
from enspara import ra
from itertools import repeat
from pathlib import Path
import matplotlib.pyplot as plt
from binding_helper import ks_weighted
import pickle


# determine weights for each sample-set taken from a given bin
# takes an array of bin weights and an array of lengths (as from RaggedArray.lengths)
# returns RaggedArray of eq_prob per bin divided by number of samples drawn from that bin.
def expand_bin_weights(eq_probs, lengths):
    return ra.RaggedArray(
        [p for i, length in enumerate(lengths) for p in repeat(eq_probs[i]/length, length)],
        lengths=lengths
    )


bins = np.arange(15, 30, 0.5)
xticks = np.arange(15, 30, 1)
fc = 'xkcd:british racing green'
obs_file = 'ntd-ctd.npy'
with open('all_ligands_experimental.pkl', 'rb') as f:
    expt = pickle.load(f)

model_type = 'mle'
tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
runtag = '12xsmina'
figdir = Path('ntd-ctd-dists')
if not figdir.is_dir():
    figdir.mkdir()
binding_name = Path(
    f'binding/resect-{resect}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
calx_name = binding_name / \
    f'extracted_scores/{runtag}/binding-calx-{model_type}'
model = Path(
    f'models/bb-all-chis-pocket-resect-{resect}/tica-{tica_lag}-msm-{msm_lag}/{k}-clusters/{model_type}/eq-probs.npy')


for mol in ['benzene', 'indole']:
    fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(12, 9))
    for rep in range(1, 4):
        rep_pre = Path(f't4l-{rep}')
        axc = axs[:, rep-1]
        axc[0].set_xticks(xticks)        
        eq_probs = np.load(str(rep_pre/model))
        reweights = ra.load(rep_pre/calx_name/f'{mol}-eq_probs.h5')
        observables = np.load(str(rep_pre/binding_name/f'receptor-observables/{obs_file}'))

        ex_eq_probs = expand_bin_weights(eq_probs, reweights.lengths)
        ks_dist, ks_prob = ks_weighted(observables, observables, ex_eq_probs.flatten(), reweights.flatten())
        axc[0].hist(observables, weights=ex_eq_probs.flatten(), bins=bins, facecolor=fc)
        axc[0].set_title('original probs')
        axc[1].hist(observables, weights=reweights.flatten(), bins=bins, facecolor=fc)
        axc[1].set_title(f'reweight by {mol}\nTwo-sided KS: diff={ks_dist:.2f}, prob={ks_prob:.4E}')
        axc[1].set_xlabel(f'Interdomain CA centroid dists ($\AA$)')
    figpath = figdir/f'{mol}-{obs_file}'
    fig.tight_layout()
    fig.savefig(str(figpath.with_suffix('.svg')), transparent=True)
    fig.savefig(str(figpath.with_suffix('.pdf')), transparent=True)
plt.show()

plt.close()