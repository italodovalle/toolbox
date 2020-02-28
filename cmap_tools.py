#! /usr/bin/python

import pandas as pd
import numpy as np
from collections import defaultdict


def parse_gct(infile):

    dt = pd.read_table(infile, sep = '\t', skiprows=2)
    cols = list(dt.columns[9:])
    desc = dt.loc[:7, cols]
    cols.append('id')
    exp_matrix = dt[cols]
    exp_matrix = dt.iloc[8:]
    attributes = defaultdict(dict)
    desc = dt.loc[:7, cols]
    for i in desc.columns[:-1]:
        for j in desc.index:
            attributes[i][desc['id'].loc[j]] = desc[i].loc[j]
    attributes = pd.DataFrame.from_dict(attributes, orient='index').reset_index()
    attributes['distil_cc_q75'] = attributes['distil_cc_q75'].astype(float)
    return(exp_matrix, attributes)
