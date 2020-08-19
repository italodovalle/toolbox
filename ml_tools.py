
import pandas as pd
import numpy as np
from collections import defaultdict
from progressbar import ProgressBar
import re
from multiprocessing import Pool
import time
import random
import os
import scipy.stats as stats
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import ShuffleSplit

def calculcate_auc(dw, vars,label = 'therapeutic',
                   ascending = True,
                   bootstrap = True, n_bootstrap = 2000,nsample = 150,
                   random_seed = 42):


    """
    Function receives a dataframe with scores and labels
    Returns auc values with confidence interval considering each score
    args:
    - vars: col names. for each col one AUC
    """

    if ascending:
        ## top predictions are the most negative values
        dw[label] = 1 - dw[label]


    res = defaultdict(dict)
    x = 0 ## counter
    for col in vars:
        sub = dw[[label, col]]
        fpr, tpr, thresholds = metrics.roc_curve(sub[label], sub[col])
        roc_auc = metrics.auc(fpr, tpr)
        res[x]['value'] = roc_auc
        res[x]['measure'] = col
        ### bootstrap
        if bootstrap:
            sub = sub.reset_index()
            rng = np.random.RandomState(random_seed)
            bootstraps = []
            for j in range(n_bootstrap):
                # bootstrap by sampling with replacement on the prediction indices
                indices = rng.random_integers(0, len(sub.index) - 1, nsample)
                boot = sub.loc[indices]
                while boot[boot[label] == 1].shape[0] == 0:
                    indices = rng.random_integers(0, len(sub.index) - 1, nsample)
                    boot = sub.loc[indices]
                fpr, tpr, thresholds = metrics.roc_curve(boot[label], boot[col])
                roc_auc_b = metrics.auc(fpr, tpr)
                bootstraps.append(roc_auc_b)
            bootstraps.sort()
            s_lower = np.percentile(bootstraps, 2.5)
            s_upper = np.percentile(bootstraps, 97.5)
            res[x]['ci_upper'] = s_upper
            res[x]['ci_lower'] = s_lower
            res[x]['error_l'] = roc_auc - s_lower
            res[x]['error_u'] = s_upper - roc_auc
        x = x + 1


    table = pd.DataFrame.from_dict(res,orient='index')
    table['chemical'] = chemical
    return (table)
