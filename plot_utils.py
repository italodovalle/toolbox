

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
from collections import defaultdict
from progressbar import ProgressBar
import scipy.stats as stats

def get_grouped_barchart(dm, id_var, value_var_1, value_var_2,
                         std_var_1, std_var_2, label_var_1, label_var_2,
                         width=0.35, n_groups = 2):

    """
    Return grouped (2 categories) barchart with error bars

    Parameters:
    -----------
    dm: pandas dataframe
        Pandas dataframe containing values to be plotted in different columns
        It also contains a column with the standard deviation
    id_var: str
        Variable that should be in the x axis
    value_var_#: str
        Column name for the element # of barchart group
    std_var_#: str
        Column name for std of the element # of the barchart group
    label_var_#: str
    width = float
    n_groups = int

    """



    fig, ax = plt.subplots()
    N = len(set(dm[id_var]))

    ind = np.arange(N)  # the x locations for the groups
           # the width of the bars

    ##groups
    ngroups = n_groups

    buf = 0
    rects = []

    groups = [[value_var_1, std_var_1, label_var_1, colors[0]],
              [value_var_2, std_var_2, label_var_2, colors[1]]]


    xticks = list(dm[id_var])
    for group in groups:

        value_var = group[0]
        std_var = group[1]
        label = group[2]
        col = group[3]

        avgs = dm[value_var].values
        avgs = avgs.reshape(1,-1).flatten()
        error_upper = (dm[value_var].values + dm[std_var].values) - dm[value_var].values
        error_lower = dm[value_var].values - (dm[value_var].values - dm[std_var].values)
        error = [error_lower.reshape(1,-1).flatten(),error_upper.reshape(1,-1).flatten()]
        rect = ax.bar(ind + buf, avgs, width, color=col, yerr=error,label = label)
        buf = buf + width

    ax.set_xticks(ind + width / 2);
    ax.set_xticklabels(xticks);
    ax.legend(loc='best')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90);
    return (ax)



def get_barplot_error (dm, id_var, value_var, std_var_1, std_var_2, pval,
                       **kwargs):


    if 'pval' in kwargs.keys():
        tr_p = kwargs['pval']
    else:
        tr_p = 0.05

    fig, ax = plt.subplots()
    N = len(set(dm[id_var]))
    ind = np.arange(N)
    width=0.35


    buf = 0
    rects = []

    xticks = list(dm[id_var])
    pvalues = list(dm[pval])

    avgs = dm[value_var].values
    avgs = avgs.reshape(1,-1).flatten()
    error_upper = np.asarray(dm[std_var_2].values - dm[value_var].values)
    error_lower = np.asarray(dm[value_var].values - dm[std_var_1].values)
    error = [error_lower.reshape(1,-1).flatten(),error_upper.reshape(1,-1).flatten()]
    rects = ax.bar(ind + buf, avgs, width,yerr=error)
    buf = buf + width



    for i,rect in enumerate(rects):
        if pvalues[i] < 0.05:
            w = ind[i]
            yloc = rect.get_height() + error[1][i] + error[1][i]/2
            xloc = w
            clr = 'black'
            align = 'left'
            ax.annotate("*", xy=(w, yloc),xytext=(-4, 0),
                                textcoords="offset points",
                                ha=align, va='center',
                                color=clr, weight='bold', fontsize = 20,clip_on=True)

    ax.set_xticks(ind + width / 2);
    ax.set_xticklabels(xticks);
    ylim = ax.get_ylim()[1] + 0.2
    ax.set_ylim(0, ylim)
    ax.legend(loc='best')
    if 'label' in kwargs.keys():
        ax.set_ylabel(kwargs['label'])

    return (fig, ax)
    
