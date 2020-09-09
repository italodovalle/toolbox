#! /usr/bin/env python

import pandas as pd
import numpy as np

from collections import defaultdict


def load_list_file(file,index_key, sep = '\t',
                 transform = None):

    """
    load list file
    first item is key, other items are elements (different lengths)
        - index_key: after split, index of the key. Elements will be 'index_key+1:'
    """

    key2val = {}
    for i in open(file).readlines():
        v = i.rstrip().split(sep)
        key = v[index_key]
        values = v[index_key+1:]
        key2val[key] = values
        if transform:
            key2val[key] = list(map(transform, values))

    return (key2val)


def df_to_dict(file, key_col, val_col, sep = ',',
               format_dict = 'dict'):

    """
    Organize data in dataframe as dict
    args:
        - sep
        - format_values: list, normal_dict

    """

    if type(file) == str:
        df = pd.read_csv(file, sep = sep)

    if type(file) == pd.core.frame.DataFrame:
        df = file

    df = df[[key_col, val_col]].drop_duplicates()


    df = df[~df[key_col].isnull()]
    df = df[~df[val_col].isnull()]

    if format_dict == 'list':
        key2val = defaultdict(list)
        for i in df.index:
            key2val[df[key_col].loc[i]].append(df[val_col].loc[i])
    elif format_dict == 'dict':
        key2val = {}
        for i in df.index:
            key2val[df[key_col].loc[i]] = df[val_col].loc[i]


    return (key2val)



def load_drug_targets(file):

    ## delete

    drug2target = {}
    for line in open(file).readlines():
        v = line.rstrip().split('\t')
        targets = v[1:]
        targets = list(map(int, targets))
        drug2target[v[0]] = targets

    return (drug2target)


def write_excel(df, outfile='out.xlsx',
                description=None):

    """
    write an excel file containing a sheet for README
    """

    if description:

        x = {'REAME':description.split('\n')}
        x = pd.DataFrame.from_dict(x, orient='columns')
    else:
        x = pd.DataFrame()

    with pd.ExcelWriter(outfile) as writer:
        df.to_excel(writer, sheet_name='data')
        x.to_excel(writer, sheet_name='README')
