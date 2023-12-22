import numpy as np
from enspara import ra
from itertools import repeat
import matplotlib.pyplot as plt
from pathlib import Path
import json
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm

# determine weights for each sample-set taken from a given bin
# takes an array of bin weights and an array of lengths (as from RaggedArray.lengths)
# returns RaggedArray of eq_prob per bin divided by number of samples drawn from that bin.


def expand_bin_weights(eq_probs, lengths):
    return ra.RaggedArray(
        [p for i, length in enumerate(lengths)
         for p in repeat(eq_probs[i]/length, length)],
        lengths=lengths
    )

# compute reweightings.


def reweighted_frames(frame_weights, kas, conc_ligand=1):
    unnormed_weights = frame_weights * (1 + kas * conc_ligand)
    return unnormed_weights / np.sum(unnormed_weights)


def ka_from_kcal_mol(fe, rt):
    return np.exp(fe / rt)**-1


minconc = 10**-9
maxconc = 10 ** -1
nconcs =  18
concentrations = np.flip(np.geomspace(minconc, maxconc, nconcs))
yticklabels = np.flip(np.geomspace(minconc, maxconc, int(nconcs/2)))
yticks = np.arange(0, nconcs, 2) +1.5
log_prob = True
bins = np.arange(0, 12, 0.5)
xticks = np.arange(0, 12, 1)
xticklabels = xticks.copy()
xtickstride = (xticks[1]-xticks[0])/2
xticks = xticks*2- xtickstride
tica_lag = 500
msm_lag = 2000
k = 75
nframes = 20
resect = '1.0'
model_type = 'mle'
runtag = '12xsmina'
if log_prob:
    figdir = Path('ligand-spyrmsd-titrations-logprob')
else:
    figdir = Path('ligand-spyrmsd-titrations')
if not figdir.is_dir():
    figdir.mkdir()
analysis_tag = figdir.stem
binding_name = Path(
    f'binding/resect-{resect}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}')
extracts_name = binding_name / f'extracted_scores/{runtag}'
calx_name = binding_name / \
    f'extracted_scores/{runtag}/binding-calx-{model_type}'
calx_fp = ('t4l-1'/calx_name/'calx.json').open('r')
calx = json.load(calx_fp)
calx_fp.close()
rt = calx['log']['rt']

fc = 'xkcd:british racing green'
molsets = {
    'benzene': '181L',
    # 'indole': '185L',
    # '1_2-dichlorobenzene': '2OTY',
    # 'benzofuran': '182L',
    # 'ethylbenzene': '7L3H',
    # 'indene': '183L',
    # 'indole': '185L',
    # 'iodobenzene': '7L3B',
    # 'iodopentafluorobenzene': '3DN3',
    # 'isobutylbenzene': '184L',
    # 'n-butylbenzene': '186L',
    # 'N-methylaniline': '2OTZ',
    # 'o-xylene': '7L3C',
    # 'propylbenzene': '7L3I',
    # 'p-xylene': '187L',
    'toluene': '7L39'
}
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
model = Path(
    f'models/bb-all-chis-pocket-resect-{resect}/tica-{tica_lag}-msm-{msm_lag}/{k}-clusters/{model_type}/eq-probs.npy')

for mol, pdb in molsets.items():
    # Set up figure and image grid
    fig = plt.figure(figsize=(9.75, 3))

    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                     nrows_ncols=(1, 3),
                     axes_pad=0.15,
                     share_all=True,
                     cbar_location="right",
                     cbar_mode="single",
                     cbar_size="7%",
                     cbar_pad=0.15)

    for rep in range(1, 4):
        rep_pre = Path(f't4l-{rep}')
        calx_fp = (rep_pre/calx_name/'calx.json').open('r')
        calx = json.load(calx_fp)
        calx_fp.close()
        rt = calx['log']['rt']
        macro_kd = calx['results'][mol]['msm K_D']**-1
        # for higher Y == lower conc, use this formula
        # ypos_kd = (np.log10(macro_kd) - np.log10(minconc))/(np.log10(maxconc)-np.log10(minconc))*nconcs
        # for higher Y == higher concentration, use this formula
        ypos_kd = (np.log10(maxconc) - np.log10(macro_kd) )/(np.log10(maxconc)-np.log10(minconc))*nconcs
        ax = grid[rep-1]
        ax.set_xlabel(f'Heavy-atom RMSD to {pdb} ligand $(\AA)$')
        ax.set_ylabel(f'[{mol}] M')
        ax.set_title(f'Replica {rep}')
        ax.plot([xticks[0], xticks[-1]+2*xtickstride], [ypos_kd, ypos_kd], label='$K_D$', color='orange')
        eq_probs = np.load(str(rep_pre/model))
        scores = ra.load(rep_pre/extracts_name/f'{mol}.h5')
        kas = ka_from_kcal_mol(scores.flatten(), rt)
        ex_eq_probs = expand_bin_weights(eq_probs, scores.lengths)
        rmsd = np.load(
            str(rep_pre/binding_name/f'noreorder-ligand-spyrmsds/{mol}/{pdb}.npy'))
        probs_vs_conc = np.zeros((len(concentrations), len(bins)-1))
        for i, conc in enumerate(concentrations):
            reweights = reweighted_frames(
                ex_eq_probs, kas, conc_ligand=conc).flatten()
            hist, _ = np.histogram(rmsd, bins=bins, weights=reweights)
            probs_vs_conc[i] = hist
        if log_prob:
            im = ax.imshow(probs_vs_conc,norm=LogNorm(vmin=1e-2, vmax=1.0))
        else:
            im = ax.imshow(probs_vs_conc, vmax=1.0, vmin=0.0)


    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.cax.colorbar(im)
    ax.cax.set_ylabel('Probability', rotation=-90, va='bottom')
    figpath = figdir/f'{mol}-{pdb}-{analysis_tag}'
    fig.tight_layout()
    fig.savefig(str(figpath.with_suffix('.pdf')), transparent=True)
    fig.savefig(str(figpath.with_suffix('.svg')), transparent=True)
plt.show()
