#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 20:54:47 2020

@author: italodovalle
"""

from . import stats_utils
from .guney_code import network_utilities
from .guney_code import wrappers
from toolbox import network_utils
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
from collections import defaultdict




class pathways:

    def __init__(self, dataframe, col_id, col_gene,
                 col_desc = None, universe=None, G=None,
                 disease_seeds=None):


        if universe:
            dataframe = dataframe[dataframe[col_gene].isin(universe)]
        pathway2entrez = defaultdict(list)
        for i in dataframe.index:
            pathway2entrez[dataframe[col_id].loc[i]].append(dataframe[col_gene].loc[i])
        self.pathway2entrez = pathway2entrez
        if col_desc:
            self.pathwayid2name = {i:j for i,j in zip(dataframe[col_id], dataframe[col_desc])}

        if G:
            self.G = G
            all_genes = list(G.nodes())
            self.all_genes = all_genes
        if disease_seeds:
            self.seeds = list(set(disease_seeds) & set(self.all_genes))
            self.disease_module = network_utils.get_lcc(self.G, self.seeds)
        #if G and disease_seeds:




    def disease_module_enrich(self,n_random = 100,
                             min_bin_size = 100, lengths = None,
                             seed=452456, ncpus = 5):

        """
        Test if pathway form a significant submodule inside the disease module
        """

        l_list  = []
        bins = network_utilities.get_degree_binning(self.G, min_bin_size, lengths)
        nodes_random = wrappers.get_random_nodes(self.seeds, self.G,
                                        bins = bins, n_random = n_random,
                                        min_bin_size = min_bin_size,
                                        seed = seed)


        all_pathways = list(self.pathway2entrez.keys())
        p = Pool(ncpus)
        #res = p.map(self.process_pathway_disease_module, all_pathways)
        res = p.map(partial(self.process_pathway_disease_module, nodes_random=nodes_random),
                    all_pathways)
        p.close()


        dic = {}
        for x in res:
            dic.update(x)

        output = pd.DataFrame.from_dict(dic, orient='index')


        return (output)



    def process_pathway_disease_module(self, pathway, nodes_random):

        dic = defaultdict(dict)

        s = list(set(self.pathway2entrez[pathway]) & set(self.all_genes))
        cc = network_utils.get_lcc(self.disease_module, s)
        real = cc.number_of_nodes()

        if real > 1:

            null = []
            for node_random in nodes_random:
                random_module = network_utils.get_lcc(self.G, node_random)
                cc_random = network_utils.get_lcc(random_module, s)
                null.append(cc_random.number_of_nodes())

            null = np.array(null)
            dic[pathway]['pathway'] = pathway
            dic[pathway]['size'] = len(s)
            dic[pathway]['mapped'] = real
            dic[pathway]['p_emp'] = null[null > real].shape[0]/null.shape[0]
            dic[pathway]['avg'] = null.mean()
            dic[pathway]['std'] = null.std()
            if null.std() == 0:
                dic[pathway]['z'] = float('nan')
            else:
                dic[pathway]['z'] = (real - null.mean())/null.std()


        return (dic)



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
