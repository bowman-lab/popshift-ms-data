import numpy as np
from enspara import ra
from itertools import repeat
import matplotlib.pyplot as plt
from pathlib import Path


# determine weights for each sample-set taken from a given bin
# takes an array of bin weights and an array of lengths (as from RaggedArray.lengths)
# returns RaggedArray of eq_prob per bin divided by number of samples drawn from that bin.
def expand_bin_weights(eq_probs, lengths):
    return ra.RaggedArray(
        [p for i, length in enumerate(lengths) for p in repeat(eq_probs[i]/length, length)],
        lengths=lengths
    )


bins = np.arange(0, 12, 0.5)
xticks = np.arange(0, 12, 1)
tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
model_type = 'mle'
runtag = '12xsmina'
figdir = Path('ligand-spyrmsd-histograms')
if not figdir.is_dir():
    figdir.mkdir()
analysis_tag = figdir.stem
binding_name = Path(f'binding/resect-{resect}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
calx_name = binding_name / f'extracted_scores/{runtag}/binding-calx-{model_type}'


fc = 'xkcd:british racing green'
# pdb = '185L'
# molsets = {
#     'benzene': '181L', 
#     'indole': '185L',
#     '1_2-dichlorobenzene': '2OTY',
#     'benzofuran': '182L',
#     'ethylbenzene': '7L3H',
#     'indene': '183L',
#     'indole': '185L',
#     'iodobenzene': '7L3B',
#     'iodopentafluorobenzene': '3DN3',
#     'isobutylbenzene': '184L',
#     'n-butylbenzene': '186L',
#     'N-methylaniline': '2OTZ',
#     'o-xylene': '7L3C',
#     'propylbenzene': '7L3I',
#     'p-xylene': '187L',
#     'toluene': '7L39'
# }
molsets = {
    'benzene': '181L', 
    'toluene': '7L39',
    'n-butylbenzene': '186L'
}
model = Path(f'models/bb-all-chis-pocket-resect-{resect}/tica-{tica_lag}-msm-{msm_lag}/{k}-clusters/{model_type}/eq-probs.npy')

for mol, pdb in molsets.items():
    fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(8,6))
    for rep in range(1, 4):
        rep_pre = Path(f't4l-{rep}')
        axc = axs[:, rep-1]
        axc[0].set_xticks(xticks)
        eq_probs = np.load(str(rep_pre/model))
        reweights = ra.load(rep_pre/calx_name/f'{mol}-eq_probs.h5')
        rmsd = np.load(str(rep_pre/binding_name/f'ligand-spyrmsds/{mol}/{pdb}.npy'))

        ex_eq_probs = expand_bin_weights(eq_probs, reweights.lengths)
        axc[0].hist(rmsd, weights=ex_eq_probs.flatten(), bins=bins, facecolor=fc)
        axc[0].set_title('original probs')
        axc[1].hist(rmsd, weights=reweights.flatten(), bins=bins, facecolor=fc)
        axc[1].set_title(f'reweight by {mol}')
        axc[1].set_xlabel(f'Heavy-atom RMSD to {pdb} ligand')
    figpath = figdir/f'{mol}-{pdb}-{analysis_tag}'
    fig.tight_layout()
    fig.savefig(str(figpath.with_suffix('.pdf')), transparent=True)
    fig.savefig(str(figpath.with_suffix('.svg')), transparent=True)
plt.show()
