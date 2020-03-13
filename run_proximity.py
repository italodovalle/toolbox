#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import random
import os
import time
from functools import partial
import re

from toolbox.guney_code import wrappers
from toolbox.guney_code import network_utilities
import toolbox.network_utils as network_utils

import logging

### Interactome params


infolder = '/home/italodovalle/flavonoids/data/'
infile = infolder + 'hiunion_interactions.csv'
header=True
sep = ','
lcc = True
columns = ['proteinA', 'proteinB']

outfolder = '/home/italodovalle/runs/polyphenols/hiunion/'
disease_genes_file = infolder + 'Guney2016_GenesDisease.tsv'
polyphenl_targets_file = infolder + 'PhenolExplorer_CTD-Stitch.csv'


ncpus = 15
n_random = 1000

outdir_tmp_files = outfolder + 'tmp'
final_outfile = outfolder + 'hiunion_zscore_proximity.csv'


def get_zscores (disease_chemical, n_random, outdir,
                     G, chemical2genes, disease2genes):

        disease, chemical = disease_chemical
        chemical = chemical.lower()
        nodes_from = set(chemical2genes[chemical]) & set(G.nodes())
        nodes_to = set(disease2genes[disease]) & set(G.nodes())
        min_bin_size = 2 * max(len(nodes_from), len(nodes_to))

        dic = network_utils.calculate_proximity_italo(G, nodes_from,nodes_to,
                                                      n_random = n_random)
        table = {}
        table['disease'] = disease
        table['n_mapped_disease'] = len(set(nodes_to))
        table['n_mapped_chemical'] = len(set(nodes_from))
        table['chemical'] = chemical
        if dic:
            table['shortest'] = dic['shortest']
            table['closest'] = dic['closest']
            table['z_shortest'] = dic['z_shortest']
            table['z_closest'] = dic['z_closest']
        else:
            table['shortest'] = float('nan')
            table['closest'] = float('nan')
            table['z_shortest'] = float('nan')
            table['z_closest'] = float('nan')

        df = pd.DataFrame.from_dict(table, orient='index').T

        if outdir:
            out_disease = re.sub('[^A-Za-z0-9]+', '', disease)
            out_chemical = re.sub('[^A-Za-z0-9]+', '', chemical)
            outname = '%s_%s'%(out_disease, out_chemical)
            outdir = os.path.abspath(outdir)
            df.to_csv(outdir + '/%s.csv'%outname)

        return(df)

def run_proximity (G, disease2genes, chemical2genes, ncpus = 15,
                   n_random = 10, outdir=None, test_run = False,
                   sp = None, node2index = None):


    finished = []
    ## retrive pairs that had their calculations already done
    #if outdir:
    #    fs = []
    #    for i in os.listdir(outdir_tmp_files):
    #        if i.endswith('.csv'):
    #            x = pd.read_csv(outdir_tmp_files + '/' + i, index_col = 0)
    #            fs.append(x)

    #    if len(fs) > 0:
    #        ds = pd.concat(fs)
    #        finished = [(i,j) for i,j in zip(ds.disease, ds.chemical)]

    samples = []
    for disease in disease2genes.keys():
        for molecule in chemical2genes.keys():
            pair = (disease,molecule)
            if not pair in finished:
                samples.append((disease,molecule))
    print ('%d Chemical-disease pairs'%len(samples))


    if test_run:
        samples = samples[:ncpus]

    p = Pool(ncpus)
    res = p.map(partial(get_zscores,
                  n_random = n_random, outdir=outdir,
                  G = G, chemical2genes = chemical2genes,
                  disease2genes = disease2genes), samples)
    p.close()

    df = pd.concat(res)
    return(df)

if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',level=logging.INFO)


    G = network_utils.parse_interactome(infile, sep, header, columns, lcc=lcc)



    logging.info('loading network: %d nodes, %d edges'%(len(G.nodes()),
                                                        len(G.edges())))

    ## create dictionary
    ## disease names are keys and elements are disease genes

    disease2genes = {}
    for i in open(disease_genes_file).readlines():
        v = i.rstrip().split('\t')
        disease = v[1]
        genes = v[2:]
        if len(genes) > 19:
            disease2genes[disease] = [int(i) for i in genes]

    ## create dictionary
    ## polyphenol names are keys and elements are polyphenol targets


    polyphenol = pd.read_csv(polyphenl_targets_file,index_col = 0)
    #polyphenol = polyphenol[(polyphenol.experimental > 0) | (polyphenol.database > 0)]
    polyphenol = polyphenol[(polyphenol.experimental > 0)]
    chemical2genes = defaultdict(list)
    for i in polyphenol.index:
        name = polyphenol.chemical.loc[i]
        chemical2genes[name].append(polyphenol.entrezid.loc[i])

    for i in chemical2genes.keys():
        chemical2genes[i] = list(set(chemical2genes[i]))



    logging.info('Test Run')

    s = time.time()
    ncpus = 10
    df = run_proximity(G, disease2genes, chemical2genes, ncpus = ncpus,
                       n_random = n_random, outdir=outdir_tmp_files, test_run=True)
    e = time.time() - s


    estimated_time = (e * (len(disease2genes) * len(chemical2genes)))/ncpus/3600

    logging.info('estimated end time: %f hours'%estimated_time)


    logging.info('Running analysis')


    s = time.time()
    df = run_proximity(G, disease2genes, chemical2genes, ncpus = ncpus,
                       n_random = n_random, outdir=outdir_tmp_files,
                       test_run=False)
    e = time.time() - s
    e = e/s

    logging.info('Finished: %f hours'%e)

    df.to_csv(final_outfile)
