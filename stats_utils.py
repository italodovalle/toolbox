import pandas as pdb
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def hypergeom_test(success_in_sample, universe, success_in_universe, sample_size):
    """success_in_sample, universe, success_in_universe, sample"""
    return (1 - stats.hypergeom.cdf(success_in_sample - 1, universe,
            success_in_universe, sample_size))


def pval_adjustment(pvals, alpha=0.05, method = 'fdr_bh'):
    adj_pval = multipletests(pvals, alpha=0.05, method=method)[1]
    return(adj_pval)

def normalized_by_rowsum(m):
    return(m.div(m.sum(axis=1), axis=0))
