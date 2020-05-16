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
                        organism = 'Homo sapiens'):
    
    
    r = pd.read_csv(
        pkg_resources.resource_filename(__name__, 
        'data/reactome/NCBI2Reactome_All_Levels.txt'),
                    sep = '\t', header=None)
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




def convert_uniprot_entrez(geneset, from_type = 'entrez'):
    
    r = pd.read_csv(
        pkg_resources.resource_filename(__name__, 
        'data/gene_ids/uniprot/HUMAN_9606_idmapping.dat'),
                    sep = '\t', header=None)
    
    if from_type == 'entrez':
        
        
    
    
    
    
    
    
    
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