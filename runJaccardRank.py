#! /usr/bin/python

import pandas as pd
import numpy as np

from collections import defaultdict
from progressbar import ProgressBar
import time

import networkx as nx
import random

import os
from scipy import sparse
from multiprocessing import Pool

import toolbox.network_utils as network_utils
import toolbox.databases_utils as databases_utils

infile = 'HI_ChengNC2018.txt'
dt = pd.read_table(infile, sep = '\t')
dt.columns = ['source', 'target', 'db']
dt.source = [int(i) for i in dt['source']]
dt.target = [int(i) for i in dt['target']]
edges = [(i,j) for i,j in zip(dt['source'],dt['target'])]

g = nx.Graph()
g.add_edges_from(edges)
G = list(nx.connected_component_subgraphs(g))[0]
print (len(G.nodes()), len(G.edges()))


seeds = pd.read_csv('seed_list.csv',index_col = 0)
seeds = list(set(seeds.entrez))
seeds = list(set(seeds) & set(G.nodes()))

def get_matrix_jc(G):
    ## get all pairwise jaccard indexes
    nodes = list(G.nodes())
    adj = nx.adj_matrix(G)
    shared_nei = adj * adj
    d = shared_nei.diagonal()

    ## matrix i,j sum degree i,j
    ## jc = overlap / (i) + (j) - overlap
    ## the steps bellow get the denominator
    a = np.tile(d,(d.shape[0],1))
    a = a.T
    b = np.tile(d,(d.shape[0],1))
    sumk = a + b
    t = shared_nei.todense()
    denom = sumk - shared_nei

    ## jaccard
    jc = t/denom

    df = pd.DataFrame(jc)
    df.columns = nodes
    df.index = nodes

    return(df)


jc = get_matrix_jc(G)


randomset = network_utils.get_random_degree_matching(G, seeds,min_bin_size=200,
n_random=1000)



def get_null(rset):
    null = {}
    for i in jc.index:
        null[i] = jc.loc[i,rset].mean()
    return(null)


dic = {}
pbar = ProgressBar()
for i in pbar(jc.index):
    dic[i] = jc.loc[i, seeds].mean()


p = Pool(15)
res = p.map(get_null,randomset)
p.close()

nullref = {i:j for i,j in zip(range(len(res)),res)}

nulldf = pd.DataFrame.from_dict(nullref, orient='columns')
print (nulldf.shape)


nulldf = pd.concat([nulldf.mean(axis=1),nulldf.std(axis=1)],axis=1)
nulldf.columns = ['avg', 'sd']


df = pd.DataFrame.from_dict(dic, orient='index')
df.columns = ['jc']
merged = pd.merge(df, nulldf, left_index=True, right_index=True)


merged['z-score'] = (merged['jc'] - merged['avg'])/merged['sd']


merged.to_csv('jaccard_index_ranking.csv')
