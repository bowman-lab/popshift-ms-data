import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv
from scipy.stats import sem

from binding_helper import plot_scatter, compare_extrema

markers = [
    'o',
    'v',
    's',
    'p'
]
pointsize = 10

with open('all_ligands_experimental.pkl', 'rb') as f:
    expt = pickle.load(f)
expt_names_lower = {i.lower(): i for i in expt.keys()}
abfe_csv = 'T4-Lysozyme_ABFE - 2022-12-23.csv'
widefig = 10
ncols = 2
nrows = 2
nframes = 20 
k = 75
tica_lag = 500
msm_lag = 2000
xtal = '186L'
resect_tag = 'resect-1.0'
model_tag = 'mle'
runtag = "12xsmina"
calx_subp = Path(f'binding/{resect_tag}/tica-{tica_lag}-msm-{msm_lag}-k-{k}-nframes-{nframes}/extracted_scores/{runtag}/calx-tica-{tica_lag}-msm-{msm_lag}-{model_tag}/calx.json')

ignore_calxtypes = {'popshift K_D', 'simple avg', 'weighted avg', 'popshift best per cluster'}

cherrypick = {'iodopentafluorobenzene', 'iodobenzene', 'N-methylaniline', '1_2-dichlorobenzene', '1-methylpyrrole'}
cherrypicktag = 'biochem_95'
filetag = '3reps-smina-4panels'
cherries = [i for i in expt if i not in cherrypick]
expt_full = np.array([expt[lig] for lig in expt])
expt_cherrypick = np.array([expt[lig] for lig in cherries])
abfe_rows = {}
with open(abfe_csv, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        compound_name = row['Compound'].replace('*','').replace(',', '_')
        if compound_name in expt:
            abfe_rows[compound_name] = row
        elif compound_name in expt_names_lower:
            abfe_rows[expt_names_lower[compound_name]] = row
        else:
            print(compound_name, f'not found in {abfe_csv}.')

for cherrypicking in [True, False]:
    fig, axs = plt.subplots(nrows, ncols, figsize=(widefig, widefig+1), sharey=True, sharex=True)
    total_x_extrem, total_y_extrem =np.array([0,-10]), np.array([0, -10])  
    reps = {}
    if cherrypicking:
        xa = expt_cherrypick
        ligs = cherries
    else:
        xa = expt_full
        ligs = list(expt.keys())
    
    # extract and average the four simulated datasets
    for rep in range(1,4):
        repp = Path(f't4l-{rep}')
        with (repp/calx_subp).open() as f:
            reps[str(repp)] = json.load(f)
    calxtypes = [key for key in reps[str(repp)]['results'][next(iter(expt))].keys() if key not in ignore_calxtypes]
    for calxind, calxt in enumerate(calxtypes):
        try:
            ax = axs.flatten()[calxind]
        except TypeError:
            ax = axs
        redock_extracts = []
        for repn, calx in reps.items():
            # rt = calx['log']['rt']
            ys = []
            for lig in expt:
                if lig in cherrypick and cherrypicking:
                    continue
                else:
                    # exptval = kcal_mol_from_kd(expt[lig] * 1e-6, rt)
                    print(repn, lig, calxt)
                    simval = calx['results'][lig][calxt]
                    ys.append(simval)
            redock_extracts.append(np.array(ys))



        ya = np.array(redock_extracts)
        ym = np.mean(ya, axis=0)
        yerr = sem(ya, axis=0)
        
        xlabel = 'Experimental $\mathrm{\Delta G_{bind}^0}$ (kcal/mol)'
        # if i == axs.shape[0] - 1:
        #     xlabel = 'Experimental $\mathrm{K_D}$ (kcal/mol)'
        # else:
        #
        #     xlabel = None
        if calxind == 0:
            ylabel = 'Docked $\mathrm{\Delta G_{bind}^0}$ (kcal/mol)'
        else:
            ylabel = None

        x_extrem, y_extrem = plot_scatter(xa, ym, ax, yerr=yerr, include_corr=True, include_rho=True, include_rmse=True,
                                          title=f'{calxt}: mean of 3 sim-sets\n', xlabel=xlabel, ylabel=ylabel, 
                                          pointshape=markers[calxind], pointsize=pointsize)
        total_x_extrem = compare_extrema(x_extrem, total_x_extrem)
        total_y_extrem = compare_extrema(y_extrem, total_y_extrem)
        # if cherrypicker:
        #     fig.savefig('correlate-expt-{}.pdf'.format('_'.join(calxt.split())), transparent=True)
        #     fig.savefig('correlate-expt-{}.svg'.format('_'.join(calxt.split())), transparent=True)
        # else:
        #     fig.savefig('cherrypicked-corr-expt-{}.svg'.format('_'.join(calxt.split())), transparent=True)
        #     fig.savefig('cherrypicked-corr-expt-{}.pdf'.format('_'.join(calxt.split())), transparent=True)


    
    calxp = Path(f'xtals/{xtal}-docks/extracts.json')
    try:
        ax = axs.flatten()[calxind+1]
    except TypeError:
        pass
    with calxp.open() as f:
        calx = json.load(f)

    ys = []
    abfes = []
    abfe_errors = []
    abfe_expt_pairs = []
    for lig in ligs:
        exptval = expt[lig]
        simval = calx[lig]
        abfe = abfe_rows[lig]['ΔGcalc corrected (kcal/mol)']
        abfe_err = abfe_rows[lig]['ΔGcalc corrected Error (kcal/mol)']
        try:
            abfes.append(float(abfe))
            abfe_errors.append(float(abfe_err))
            abfe_expt_pairs.append(exptval)
        except ValueError:
            print(lig, 'had no ABFE.')
        ys.append(simval)

    ya = np.array(ys)
    abfes = np.array(abfes)
    abfe_errors = np.array(abfe_errors)
    abfe_xs = np.array(abfe_expt_pairs)

    x_extrem, y_extrem = plot_scatter(xa, ya, ax, include_corr=True, include_rho=True, include_rmse=True,
                                    xlabel=xlabel, ylabel=ylabel, pointshape=markers[calxind+1], pointsize=pointsize,
                                title='Docking to n-butylbenzene holo structure,\n')
    try:
        ax = axs.flatten()[calxind+2]
    except TypeError:
        pass
    x_extrem, y_extrem = plot_scatter(abfe_xs, abfes, ax, yerr=abfe_errors, include_corr=True, include_rho=True, include_rmse=True,
                                    xlabel=xlabel, ylabel=ylabel, pointshape=markers[calxind+2], pointsize=pointsize,
                                title='Absolute Binding free energies (dissapearing ligand),\n')
    total_x_extrem = compare_extrema(x_extrem, total_x_extrem)
    total_y_extrem = compare_extrema(y_extrem, total_y_extrem)


    ax.set_xlim(total_x_extrem[0], total_x_extrem[1])
    ax.set_ylim(total_y_extrem[0], total_y_extrem[1])
    # ax.legend()
    # axs[0].set_xlim(total_x_extrem[0], total_x_extrem[1])
    # axs[0].set_ylim(total_y_extrem[0], total_y_extrem[1])
    # axs[0].legend()

    fig.tight_layout()
    if cherrypicking:
        fig.savefig(f'dgs_{filetag}_vs_expt_{cherrypicktag}.svg', transparent=True)
        fig.savefig(f'dgs_{filetag}_vs_expt_{cherrypicktag}.pdf', transparent=True)
    else:
        fig.savefig(f'dgs_{filetag}_vs_expt.pdf', transparent=True)
        fig.savefig(f'dgs_{filetag}_vs_expt.pdf', transparent=True)

plt.show()
