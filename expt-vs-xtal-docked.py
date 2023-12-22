import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from binding_helper import plot_scatter, compare_extrema


with open('all_ligands_experimental.pkl', 'rb') as f:
    expt = pickle.load(f)

widefig = 9
ncols = 1
nrows = 1

ligands = [i for i in expt]
xtal = '186L'
not_b_95 = {'iodopentafluorobenzene', 'iodobenzene',
            'N-mythlaniline', '1_2-dichlorobenzene', '1-methylpyrrole'}
b_95 = [i for i in ligands if i not in not_b_95]
filetag = 'biochem-95-vs-xtal-docks'
calxp = Path(f'xtals/{xtal}-docks/extracts.json')
with calxp.open() as f:
    calx = json.load(f)

fig, ax = plt.subplots(nrows, ncols, figsize=(widefig, widefig/ncols), sharey=True, sharex=True)
total_x_extrem, total_y_extrem = np.array([0,-10]), np.array([0, -10])

xs = []
ys = []
for lig in b_95:
    exptval = expt[lig]
    simval = calx[lig]
    ys.append(simval)
    xs.append(exptval)

xa = np.array(xs)
ya = np.array(ys)
ylabel = 'Docked $\mathrm{\Delta G_{bind}^0}$ (kcal/mol)'
xlabel = 'Experimental $\mathrm{\Delta G_{bind}^0}$ (kcal/mol)'

x_extrem, y_extrem = plot_scatter(xa, ya, ax, include_corr=True, include_rho=True, include_rmse=True,
                                  xlabel=xlabel, ylabel=ylabel,
                            title='Docking to n-butylbenzene holo structure')
total_x_extrem = compare_extrema(x_extrem, total_x_extrem)
total_y_extrem = compare_extrema(y_extrem, total_y_extrem)


ax.set_xlim(total_x_extrem[0], total_x_extrem[1])
ax.set_ylim(total_y_extrem[0], total_y_extrem[1])
ax.legend()

fig.tight_layout()
fig.savefig(f'{filetag}.pdf', transparent=True)
fig.savefig(f'{filetag}.pdf', transparent=True)

plt.show()
