#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 20:54:47 2020

@author: italodovalle
"""

from . import stats_utils
import pandas as pd
from collections import defaultdict
import pkg_resources


def enrichment_reactome(geneset, alpha = 0.05,correction = 'fdr_bh',
                        organism = 'Homo sapiens',
                        pathway_file=None,label=None):

    """
    Pathway Enrichment Analysis (Reactome)

    geneset: list of entrez ids
    """

    if type(pathway_file) == pd.core.frame.DataFrame:
        r = pathway_file
    elif type(pathway_file) == str:
        r = pd.read_csv(pathway_file, sep = '\t', header=None,
        low_memory=False)
    else:
        url = 'https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt'
        r = pd.read_csv(url, sep = '\t', header=None,
        low_memory=False)
    r = r[r[5] == organism]


    geneset = list(map(int, geneset))
    geneset = list(map(str, geneset))

    universe = list(set(r[0]))
    universe = list(map(str, universe))

    pathways = list(set(r[1]))

    pathway2name = {i:j for i,j in zip(r[1],r[3])}

    table = defaultdict(dict)
    c = 0

    for pathway in pathways:

        sample = list(set(r[r[1] == pathway][0]))
        success_in_universe = list(set(geneset) & set(universe))
        success_in_sample = list(set(sample) & set(geneset))
        pval = stats_utils.hypergeom_test(len(success_in_sample), len(universe),
                        len(success_in_universe), len(sample))

        table[c]['PathwayID'] = pathway
        table[c]['PathwayName'] = pathway2name[pathway]
        table[c]['gene_ratio'] = 1. * len(success_in_sample)/len(sample)
        table[c]['pvalue'] = pval
        c = c + 1

    table = pd.DataFrame.from_dict(table, orient='index')

    table['pvalue_adj'] = stats_utils.pval_adjustment(table.pvalue, alpha=alpha,
                                          method = correction)


    if label:
        table['label'] = label

    return (table)




def convert_uniprot_entrez(geneset, from_type = 'entrez'):

    r = pd.read_csv(
        pkg_resources.resource_filename(__name__,
        'data/gene_ids/uniprot/HUMAN_9606_idmapping.dat'),
                    sep = '\t', header=None)

    #if from_type == 'entrez':

    r = r[r[5] == organism]

    universe = list(set(r[0]))

    pathways = list(set(r[1]))

    pathway2name = {i:j for i,j in zip(r[1],r[2])}

    table = defaultdict(dict)
    c = 0

    for pathway in pathways:

        sample = list(set(r[r[1] == pathway][0]))
        success_in_universe = list(set(geneset) & set(universe))
        success_in_sample = list(set(sample) & set(geneset))
        pval = stats_utils.enrichment_test(len(success_in_sample), len(universe),
                        len(success_in_universe), len(sample))

        table[c]['PathwayID'] = pathway
        table[c]['PathwayName'] = pathway2name[pathway]
        table[c]['gene_ratio'] = 1. * len(success_in_sample)/len(sample)
        table[c]['pvalue'] = pval
        c = c + 1

    table = pd.DataFrame.from_dict(table, orient='index')

    table['pvalue_adj'] = stats_utils.pval_adjustment(table.pvalue, alpha=alpha,
                                          method = correction)

    return (table)
