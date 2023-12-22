import numpy as np
from scipy.stats import gaussian_kde, spearmanr, distributions, pearsonr
from pathlib import Path
import matplotlib.pyplot as plt


def ks_weighted(data1, data2, wei1, wei2, alternative='two-sided'):
    ix1 = np.argsort(data1)
    ix2 = np.argsort(data2)
    data1 = data1[ix1]
    data2 = data2[ix2]
    wei1 = wei1[ix1]
    wei2 = wei2[ix2]
    data = np.concatenate([data1, data2])
    cwei1 = np.hstack([0, np.cumsum(wei1)/sum(wei1)])
    cwei2 = np.hstack([0, np.cumsum(wei2)/sum(wei2)])
    cdf1we = cwei1[np.searchsorted(data1, data, side='right')]
    cdf2we = cwei2[np.searchsorted(data2, data, side='right')]
    d = np.max(np.abs(cdf1we - cdf2we))
    # calculate p-value
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    m, n = sorted([float(n1), float(n2)], reverse=True)
    en = m * n / (m + n)
    if alternative == 'two-sided':
        prob = distributions.kstwo.sf(d, np.round(en))
    else:
        z = np.sqrt(en) * d
        # Use Hodges' suggested approximation Eqn 5.3
        # Requires m to be the larger of (n1, n2)
        expt = -2 * z**2 - 2 * z * (m + 2*n)/np.sqrt(m*n*(m+n))/3.0
        prob = np.exp(expt)
    return d, prob


def kcal_mol_from_kd(kd, rt):
    return rt * np.log(kd)

def make_gaussian_kde(x, y):
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    return x[idx], y[idx], z[idx]


def tag_from_fname(fname, splitstr='-scores'):
    p = Path(fname)
    tag, _ = p.stem.split(splitstr, 1)
    return tag.upper()


def compare_extrema(extrema_A, extrema_B):
    return np.array([min(extrema_A[0], extrema_B[0]),
                     max(extrema_A[1], extrema_B[1])])


def plot_scatter(ref_log, deriv_log, ax: plt.axis, yerr=None, xerr=None, xlabel=None, ylabel=None, title=None, include_corr=False,
                 include_rmse=False, include_rho=True, corrfmt='Corr. {:.3f} (p={:.3f})', delta_guides=1.2, ax_pad=0.5,
                 pointsize=8, pointshape='o', datalabel=None):
    one_one_x = np.array([min(ref_log.min(), deriv_log.min())-ax_pad,
                          max(deriv_log.max(), ref_log.max())+ax_pad])
    one_one_y = one_one_x
    # this will work if the points for each frame were taken in the same order,
    # so that the rows of the logs match 1-1
    x, y = ref_log, deriv_log
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        if include_corr:
            res = pearsonr(x, y)
            corr = np.corrcoef(x, y)[0,1]
            pval = res.pvalue
            title += corrfmt.format(corr, pval)
        if include_rho:
            res = spearmanr(x, y)
            rho = res.correlation
            rpval = res.pvalue
            title += ', Rho {:.3f} (p={:.3f})'.format(rho, rpval)
        if include_rmse:
            rmse = np.sqrt(np.mean((x - y)**2))
            title += ', RMSE {:.2f}'.format(float(rmse))
        ax.set_title(title)
    elif datalabel:
        if include_corr:
            corr = np.corrcoef(x, y)
            datalabel += corrfmt.format(corr[0, 1])
        if include_rho:
            rho = spearmanr(x, y).correlation
            datalabel += ', Rho {:.3f}'.format(rho)
        if include_rmse:
            rmse = np.sqrt(np.mean((x - y)**2))
            datalabel += ', RMSE {:.2f}'.format(float(rmse))
    else:
        corr = np.corrcoef(x, y)
        ax.set_title(corrfmt.format(corr[0, 1]))
    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=pointshape, ms=pointsize, lw=pointsize/4, label=datalabel, ecolor='black', alpha=0.5)
    ax.plot(one_one_x, one_one_y, color='k')
    if delta_guides:
        ax.plot(one_one_x-delta_guides, one_one_y+delta_guides, color='b')
        ax.plot(one_one_x+delta_guides, one_one_y-delta_guides, color='b')

    ax.set_xlim(one_one_x[0], one_one_x[1])
    ax.set_ylim(one_one_y[0], one_one_y[1])

    return one_one_x, one_one_y


def plot_scatter_dens(ref_log, deriv_log, ax: plt.axis, xlabel=None, ylabel=None, title=None, include_corr=False,
                      include_rho=True, corrfmt=', Corr. {:.3f}', delta_guides=1.2, ax_pad=0.5):
    one_one_x = np.array([min(ref_log.min(), deriv_log.min())-ax_pad,
                          max(deriv_log.max(), ref_log.max())+ax_pad])
    one_one_y = one_one_x
    # this will work if the points for each frame were taken in the same order,
    # so that the rows of the logs match 1-1
    x, y, z = make_gaussian_kde(ref_log, deriv_log)
    ax.scatter(x, y, c=z, s=8)
    ax.plot(one_one_x, one_one_y, color='k')
    if delta_guides:
        ax.plot(one_one_x-delta_guides, one_one_y+delta_guides, color='b')
        ax.plot(one_one_x+delta_guides, one_one_y-delta_guides, color='b')

    ax.set_xlim(one_one_x[0], one_one_x[1])
    ax.set_ylim(one_one_y[0], one_one_y[1])
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        if include_corr:
            corr = np.corrcoef(x, y)
            title += corrfmt.format(corr[0, 1])
        if include_rho:
            rho = spearmanr(x, y).correlation
            title += ', Rho. {:.3f}'.format(rho)
        ax.set_title(title)
    else:
        corr = np.corrcoef(x, y)
        ax.set_title(corrfmt.format(corr[0, 1]))
    return one_one_x, one_one_y


def plot_scatter_dens_no_rangematch(ref_log, deriv_log, ax: plt.axis, xlabel=None, ylabel=None, title=None, include_corr=False,
                      corrfmt=', Corr. {:.3f}', delta_guides=1.2, ax_pad=0.5):
    # one_one_x = np.array([min(ref_log.min(), deriv_log.min())-ax_pad,
    #                       max(deriv_log.max(), ref_log.max())+ax_pad])
    # one_one_y = one_one_x
    # this will work if the points for each frame were taken in the same order,
    # so that the rows of the logs match 1-1
    x, y, z = make_gaussian_kde(ref_log, deriv_log)
    xs = (np.min(ref_log), np.max(ref_log))
    ax.scatter(x, y, c=z, s=8)
    linear_model = np.polyfit(ref_log, deriv_log, 1)
    linear_fn = np.poly1d(linear_model)
    ax.plot(xs, linear_fn(xs), color='k')


    # if delta_guides:
    #     ax.plot(one_one_x-delta_guides, one_one_y+delta_guides, color='b')
    #     ax.plot(one_one_x+delta_guides, one_one_y-delta_guides, color='b')

    # ax.set_xlim(one_one_x[0], one_one_x[1])
    # ax.set_ylim(one_one_y[0], one_one_y[1])
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        if include_corr:
            corr = np.corrcoef(x, y)
            title += corrfmt.format(corr[0, 1])
            title += ' line: {:.3f}x + {:.3f}'.format(*linear_model)
        ax.set_title(title)
    else:
        corr = np.corrcoef(x, y)
        ax.set_title(corrfmt.format(corr[0, 1]))

# find indices where two arrays are more different than  some threshold value
def get_disagreement(x, y, threshold=1.0):
    return np.where(np.abs(x - y) > threshold)


# find indices where two arrays are more different than some signed difference threshold value.
def get_signed_disagreement(startpoint, endpoint, threshold=1.0):
    return np.where(endpoint - startpoint > threshold)

# find indices where two arrays are more different than  some threshold value
def get_agreement(x, y, threshold=1.0):
    return np.where(np.abs(x - y) < threshold)


# find indices where two arrays are more different than some signed difference threshold value.
def get_signed_agreement(startpoint, endpoint, threshold=1.0):
    return np.where(endpoint - startpoint < threshold)


# find all docking score indices above/below some value
def get_high_scores(scores, threshold=-6.0):
    return np.where(scores > threshold)


def get_low_scores(scores, threshold=-6.0):
    return np.where(scores < threshold)
