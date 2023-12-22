import matplotlib.pyplot as plt
from scipy.stats import sem
import numpy as np


def plot_one_vampseries(ax, k, vamp_fn, color, label, fmt='o'):
    vamp_series = np.load(vamp_fn)
    vamp_means = np.nanmean(vamp_series, axis=1)
    vamp_errs = sem(vamp_series, axis=1)
    ax.errorbar(k, vamp_means, yerr=vamp_errs, color=color, label=label, fmt=fmt)
    return vamp_means, vamp_errs


k = [25, 50, 75, 100, 200, 300, 400, 500]
lags = [100, 200, 500, 1000]
reprange = range(1, 4)
# n_datasets = len(lags) + len(reprange)
n_datasets = len(reprange)
tts_colors = {'train': plt.cm.get_cmap('Purples')(np.linspace(0.5, 1, n_datasets)),
              'test': plt.cm.get_cmap('Greens')(np.linspace(0.5, 1, n_datasets))}
figsize_scaler = 5
figsize=(figsize_scaler*5, len(lags)*figsize_scaler)
fig, axs = plt.subplots(len(lags), 1, sharex=True, sharey=True, figsize=figsize)
for lag_ix, lag in enumerate(lags):
    ax = axs[lag_ix]
    means = {'train': [],
             'test': []}
    errs = {'train': [],
            'test': []}
    for tts in ('train', 'test'):
        for rep_ix, rep in enumerate(range(1,4)):
            vamp_fn = f't4l-{rep}/clustering/vamp-scan/t4l-{rep}-backbone-all-chis-' \
                      f'dihedrals-pocket-tica-lag-{lag}-tica-reduced-vamp-scan-' \
                      f'lagtime-{lag}-k-25-50-75-100-200-300-400-500-{tts}.npy'
            color = tts_colors[tts][rep_ix]
            label = f'{tts} lag {lag} {rep}'
            vamp_means, vamp_errs = plot_one_vampseries(ax, k, vamp_fn, color, label)
            means[tts].append(vamp_means)
            errs[tts].append(vamp_errs)
        avgmeans = np.nanmean(means[tts], axis=0)
        testerrs = np.array(errs[tts])
        quaderrs = np.sqrt(np.nansum(testerrs * testerrs, axis=0))
        ax.errorbar(k, avgmeans, yerr=quaderrs, color='k', lw=2, label=f'{tts}, rep-mean')
    ax.legend()
    ax.set_title(f'{lag} step lag')

ax.set_xticks(k)
ax.set_xticklabels(map(str, k))
fig.tight_layout()
fig.savefig('tica-vamp-scan.pdf', transparent=True)
# plt.show()
