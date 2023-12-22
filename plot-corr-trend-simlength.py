import matplotlib.pyplot as plt
from enspara import ra
import numpy as np
import json
import pickle
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, permutation_test, sem


def getrho(x, y):
    return spearmanr(x, y).statistic


def perm_pval(x, y, n_resamples=200):
    return permutation_test((x,), lambda p: getrho(p, y),
                            permutation_type='pairings',
                            n_resamples=n_resamples).pvalue


with open('all_ligands_experimental.pkl', 'rb') as f:
    all_ligands_expt = pickle.load(f)

ligands = [i for i in all_ligands_expt]
not_b_95 = {'iodopentafluorobenzene', 'iodobenzene',
            'N-mythlaniline', '1_2-dichlorobenzene', '1-methylpyrrole'}
b_95 = [i for i in ligands if i not in not_b_95]
runtag = '12xsmina'

reps = list(range(1, 4))
nreps = len(reps)
expt_arr = [all_ligands_expt[lig] for lig in b_95]
lag = 500
k = 75
msm_lag = 2000
# nframes = [1, 3, 5, 10, 20]
resects = [0.2, 0.4, 0.6, 0.8, 1.0]
nfr = 20
model_type = 'mle'
best_corr = []
best_corr_err = []
best_rho = []
best_rho_err = []
best_rmse = []
best_rmse_err = []
msm_corr = []
msm_corr_err = []
msm_rho = []
msm_rho_err = []
msm_rmse = []
msm_rmse_err = []
for rct in resects:
    rep_best_corr = []
    rep_best_rho = []
    rep_msm_corr = []
    rep_msm_rho = []
    rep_best_rmse = []
    rep_msm_rmse = []
    for rep in reps:
        top = Path(f't4l-{rep}/binding/resect-{rct}')
        extractp = top/f'tica-{lag}-msm-{msm_lag}-k-{k}-nframes-{nfr}' / \
            f'extracted_scores/{runtag}/binding-calx-{model_type}/calx.json'
        with extractp.open('r') as f:
            extract = json.load(f)
        msm_dGs = []
        best_score = []
        results = extract['results']
        for ligand in b_95:
            try:
                ligres = results[ligand]
            except KeyError as e:
                print(extractp)
                raise e
            msm_dGs.append(ligres['msm dG'])
            best_score.append(ligres['best score'])

        msm_dGs = np.array(msm_dGs)
        best_score = np.array(best_score)
        rep_best_rmse.append(np.sqrt(np.mean((best_score - expt_arr)**2)))
        rep_msm_rmse.append(np.sqrt(np.mean((msm_dGs - expt_arr)**2)))
        msm_pearson = pearsonr(msm_dGs, expt_arr)
        rep_msm_corr.append(msm_pearson.statistic)
        # msm_ci = msm_pearson.confidence_interval()
        msm_rho_s = getrho(msm_dGs, expt_arr)
        rep_msm_rho.append(msm_rho_s)
        # msm_rho_p = perm_pval(msm_dGs, expt_arr, n_resamples=300)
        best_pearson = pearsonr(best_score, expt_arr)
        rep_best_corr.append(best_pearson.statistic)
        best_rho_s = getrho(best_score, expt_arr)
        rep_best_rho.append(best_rho_s)
        # best_rho_p = perm_pval(best_score, expt_arr, n_resamples=300)
        # sel_model_df['msm rho pval'].append(msm_rho_p)
        # sel_model_df['best score rho pval'].append(best_rho_p)
    best_corr.append(np.mean(rep_best_corr))
    best_corr_err.append(sem(rep_best_corr))
    best_rho.append(np.mean(rep_best_rho))
    best_rho_err.append(sem(rep_best_rho))
    best_rmse.append(np.mean(rep_best_rmse))
    best_rmse_err.append(sem(rep_best_rmse))
    msm_corr.append(np.mean(rep_msm_corr))
    msm_corr_err.append(sem(rep_msm_corr))
    msm_rho.append(np.mean(rep_msm_rho))
    msm_rho_err.append(sem(rep_msm_rho))
    msm_rmse.append(np.mean(rep_msm_rmse))
    msm_rmse_err.append(sem(rep_msm_rmse))

fig, ax = plt.subplots()
xticks = resects
# ax.errorbar(resects, msm_corr, yerr=msm_corr_err, label='MSM binding corr.')
ax.errorbar(resects, msm_rho, yerr=msm_rho_err, label='PopShift')
# ax.errorbar(resects, best_corr, yerr=best_corr_err, label='Best score corr.')
ax.errorbar(resects, best_rho, yerr=best_rho_err, label='Best score')
ax.set_xticks(xticks)
ax.grid(visible=True, axis='y', alpha=0.3)
ax.legend()

ax.set_title(f'Correlation with experiment vs simulation length across {nreps} replicas')
ax.set_ylabel('Spearman Correlation')
# ax = axs[1]
# ax.errorbar(resects, msm_rmse, yerr=msm_rmse_err, label='MSM RMSE')
# ax.errorbar(resects, best_rmse, yerr=best_rmse_err, label='Best score RMSE')
# ax.set_xticks(xticks)
# ax.set_xlabel('Receptor Conf. Samples per MSM Bin')
# ax.legend()
# ax.set_title(f'RMSE with experiment vs bin samples across {nreps} full replicas')
# ax.set_ylabel('kcal/mol')
fig.savefig(f'corr-trend-simlength-{model_type}-{nfr}.svg', transparent=True)
fig.savefig(f'corr-trend-simlength-{model_type}-{nfr}.pdf', transparent=True)
plt.show()

