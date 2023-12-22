import numpy as np
from enspara import ra
from itertools import repeat
import matplotlib.pyplot as plt
from pathlib import Path
import pickle

# determine weights for each sample-set taken from a given bin
# takes an array of bin weights and an array of lengths (as from RaggedArray.lengths)
# returns RaggedArray of eq_prob per bin divided by number of samples drawn from that bin.
def expand_bin_weights(eq_probs, lengths):
    return ra.RaggedArray(
        [p for i, length in enumerate(lengths) for p in repeat(eq_probs[i]/length, length)],
        lengths=lengths
    )


bins = np.arange(-180,180, 10)
xticks = np.arange(-180,180, 30)
fc = 'xkcd:british racing green'

tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
model_type = 'mle'
runtag = '12xsmina'
figdir = Path('val111_chi1')
if not figdir.is_dir():
    figdir.mkdir()
binding_name = Path(
    f'binding/resect-{resect}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
calx_name = binding_name / \
    f'extracted_scores/{runtag}/binding-calx-{model_type}'
obs_file = 'val111_chi1.npy'
with open('all_ligands_experimental.pkl', 'rb') as f:
    expt = pickle.load(f)
model = Path(
    f'models/bb-all-chis-pocket-resect-{resect}/tica-{tica_lag}-msm-{msm_lag}/{k}-clusters/{model_type}/eq-probs.npy')

for mol in ['benzene', 'o-xylene']:
    fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(8,6))
    for rep in range(1, 4):
        rep_pre = Path(f't4l-{rep}')
        axc = axs[:, rep-1]
        axc[0].set_xticks(xticks)
        eq_probs = np.load(str(rep_pre/model))
        reweights = ra.load(rep_pre/calx_name/f'{mol}-eq_probs.h5')
        observables = np.load(str(rep_pre/binding_name/f'receptor-observables/{obs_file}'))

        ex_eq_probs = expand_bin_weights(eq_probs, reweights.lengths)
        axc[0].hist(observables, weights=ex_eq_probs.flatten(), bins=bins, facecolor=fc)
        axc[0].set_title('original probs')
        axc[1].hist(observables, weights=reweights.flatten(), bins=bins, facecolor=fc)
        axc[1].set_title(f'reweight by {mol}')
        axc[1].set_xlabel(f'Val111 $\chi_1$')
        axc[1].set_xticklabels(axc[1].get_xticklabels(), rotation=45)
    figpath = figdir/f'{mol}-{obs_file}'
    fig.tight_layout()
    fig.savefig(str(figpath.with_suffix('.svg')), transparent=True)
    fig.savefig(str(figpath.with_suffix('.pdf')), transparent=True)
    # fig.close()
plt.show()
