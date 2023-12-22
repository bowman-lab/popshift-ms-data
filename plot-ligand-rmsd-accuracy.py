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
        [p for i, length in enumerate(lengths)
         for p in repeat(eq_probs[i]/length, length)],
        lengths=lengths
    )


# def plot bars across an accumulated set, with class offsets
def plot_grouped_bars(ax, means_by_group, ylabel, title, xticks, xlabels, legend=None, barwidth=0.25, yvalpad=0, ylim=None):
    for i, (accuracy_class, mean_accum) in enumerate(means_by_group.items()):
        offset = barwidth * i
        rects = ax.bar(xpos + offset, mean_accum, barwidth, label=accuracy_class)
        ax.bar_label(rects, padding=yvalpad, fmt='%.2f')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(xticks+barwidth, xlabels)
    if legend:
        ax.legend(**legend)
    
    if ylim:
        ax.set_ylim(*ylim)



bins = np.arange(0, 12, 0.5)
xtal_cat = np.where(bins < 2)
alternative_cat = np.where((2 <= bins) & (bins < 4))
miss_cat = np.where(4 <= bins[:-1])
xticks = np.arange(0, 12, 1)
tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
model_type = 'mle'
runtag = '12xsmina'
figdir = Path('ligand-spyrmsd-accuracy')
if not figdir.is_dir():
    figdir.mkdir()
analysis_tag = figdir.stem
binding_name = Path(
    f'binding/resect-{resect}/'
    f'tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
calx_name = binding_name / \
    f'extracted_scores/{runtag}/binding-calx-{model_type}'


molsets = {
    'benzene': '181L',
    # 'indole': '185L',
    '1_2-dichlorobenzene': '2OTY',
    'benzofuran': '182L',
    'ethylbenzene': '7L3H',
    'indene': '183L',
    # 'indole': '185L',
    'iodobenzene': '7L3B',
    'iodopentafluorobenzene': '3DN3',
    'isobutylbenzene': '184L',
    'n-butylbenzene': '186L',
    # 'N-methylaniline': '2OTZ',
    'o-xylene': '7L3C',
    'propylbenzene': '7L3I',
    'p-xylene': '187L',
    'toluene': '7L39'
}
xpos = np.arange(len(molsets.keys()))
barwidth = 0.30


model = Path(
    f'models/bb-all-chis-pocket-resect-{resect}/'
    f'tica-{tica_lag}-msm-{msm_lag}/{k}-clusters/'
    f'{model_type}/eq-probs.npy')

fig, axs = plt.subplots(2,1, sharex=True, sharey=True, figsize=(14, 8))
xtal_tag = 'crystal-like'
alt_tag = 'alternative'
miss_tag = 'miss'
apo_means = {
    xtal_tag : [],
    alt_tag: [],
    miss_tag : []
}
holo_means = {
    xtal_tag : [],
    alt_tag : [],
    miss_tag : []
}
for mol, pdb in molsets.items():
    apo_accum = {
        xtal_tag : [],
        alt_tag: [],
        miss_tag: []
    }
    holo_accum = {
        xtal_tag: [],
        alt_tag: [],
        miss_tag: []
    }
    for rep in range(1, 4):
        rep_pre = Path(f't4l-{rep}')
        eq_probs = np.load(str(rep_pre/model))
        reweights = ra.load(rep_pre/calx_name/f'{mol}-eq_probs.h5')
        rmsd = np.load(
            str(rep_pre/binding_name/f'noreorder-ligand-spyrmsds/{mol}/{pdb}.npy'))

        ex_eq_probs = expand_bin_weights(eq_probs, reweights.lengths)
        apo_binprobs, _ = np.histogram(
            rmsd, weights=ex_eq_probs.flatten(), bins=bins)
        est_holo_binprobs, _ = np.histogram(
            rmsd, weights=reweights.flatten(), bins=bins)
        apo_accum[xtal_tag].append(np.sum(apo_binprobs[xtal_cat]))
        apo_accum[alt_tag].append(np.sum(apo_binprobs[alternative_cat]))
        apo_accum[miss_tag].append(np.sum(apo_binprobs[miss_cat]))
        
        holo_accum[xtal_tag].append(np.sum(est_holo_binprobs[xtal_cat]))
        holo_accum[alt_tag].append(np.sum(est_holo_binprobs[alternative_cat]))
        holo_accum[miss_tag].append(np.sum(est_holo_binprobs[miss_cat]))
    for k, v in apo_accum.items():
        apo_means[k].append(np.mean(v))
    for k, v in holo_accum.items():
        holo_means[k].append(np.mean(v))
# plot apo-weight bars
plot_grouped_bars(
    axs[0],
    apo_means,
    'fraction total prob.',
    'Similarity to crystal pose - Ligand-Free populations',
    xpos,
    molsets.keys(),
    legend=dict(ncols=len(apo_means.keys())),
    barwidth=barwidth
    # ylim=(0, 1)
)
plot_grouped_bars(
    axs[1],
    holo_means,
    'fraction total prob.',
    'Similarity to crystal pose - Ligand-Saturated populations',
    xpos,
    molsets.keys(),
    barwidth=barwidth
    # ylim=(0, 1)
)
axs[1].set_xticklabels(molsets.keys(), rotation=45, ha='right')
fig.tight_layout()
fig.savefig(str(figdir/'fraction-correct.pdf'), transparent=True)
fig.savefig(str(figdir/'fraction-correct.svg'), transparent=True)
# fig.autofmt_xdate(rotation=30)
plt.show()
