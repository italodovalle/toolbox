#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:06:16 2017

@author: italodovalle
"""



import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import random
from random import shuffle
import separation
from guney_code.wrappers import network_utilities, get_random_nodes, calculate_closest_distance, calculate_separation_proximity, calculate_proximity
from scipy import stats
import mygene
import os
from progressbar import ProgressBar

class NetworkUtilsInputError(Exception):
    pass


def calculate_proximity_italo(network, nodes_from, nodes_to, 
                        nodes_from_random=None, nodes_to_random=None, 
                        bins=None, n_random=1000, min_bin_size=100, 
                        seed=452456, sp = None, node2index = None,lengths = None):
    
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    last edit: Italo Oct 12, 2019
    """
    
    
    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
        return None # At least one of the node group not in network
    
    d = calculate_distances(network, nodes_from, nodes_to,sp, node2index)
    
    if n_random:
    
        if bins is None and (nodes_from_random is None or nodes_to_random is None):
            bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
        if nodes_from_random is None:
            nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        if nodes_to_random is None:
            nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        random_values_list = list(zip(nodes_from_random, nodes_to_random))
        #values = np.empty(len(nodes_from_random)) #n_random
        null = []
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            res = calculate_distances(network, nodes_from, nodes_to,sp, node2index)
            null.append(res)
            

        null_s = []
        null_c = []
        for i in range(len(null)):
            null_s.append(null[i]['shortest'])
            null_c.append(null[i]['closest'])

        d['avg_shortest'],d['std_shortest'] = np.mean(null_s), np.std(null_s)
        d['z_shortest'] = (d['shortest'] - d['avg_shortest'])/d['std_shortest']

        d['avg_closest'],d['std_closest'] = np.mean(null_c), np.std(null_c)
        d['z_closest'] = (d['closest'] - d['avg_closest'])/d['std_closest']
    
    return (d)



def distance2component(C,t,G):
    dic = defaultdict(dict)
    for s in C:
        if nx.has_path(G, s, t):
            dic[t][s] = nx.shortest_path_length(G,s,t)
        else:
            dic[t][s] = float('nan')
    d = pd.DataFrame.from_dict(dic)
    v = stats.gmean(list(d.T.iloc[0]))
    return(v)

def calculate_shortest_components (network, nodes_from, nodes_to, 
                        nodes_from_random=None, nodes_to_random=None, 
                        bins=None, n_random=1000, min_bin_size=100, 
                        seed=452456, lengths=None):
    
    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
        return None # At least one of the node group not in network
    
    d = shortest_components_geom(nodes_from, nodes_to, network)
    
    if n_random:
    
        if bins is None and (nodes_from_random is None or nodes_to_random is None):
            bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
        if nodes_from_random is None:
            nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        if nodes_to_random is None:
            nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        random_values_list = list(zip(nodes_from_random, nodes_to_random))
        #values = np.empty(len(nodes_from_random)) #n_random
        null = []
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            res = shortest_components_geom(nodes_from, nodes_to, network)
            null.append(res)
            #values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)

        #pval = float(sum(values <= d)) / len(values) # needs high number of n_random


        dic = {}
        dic['shortest'] = d 
        dic['mean'] = np.mean(null)
        dic['std'] = np.std(null)
        dic['z'] = (d - dic['mean'])/dic['std']
    
    return (d)


def parse_interactome(infile, sep='\t', header=False, columns=[], lcc = False):
    
    """
    
    """
    
    if header:
        dt = pd.read_table(infile,sep = sep)
        edges = zip(dt[columns[0]], dt[columns[1]])
    else:
        header = None
        dt = pd.read_table(infile,sep = sep,header=header)
        edges = zip(dt[0], dt[1])

    G = nx.Graph()
    G.add_edges_from(edges)

    if lcc:
        g = list(nx.connected_component_subgraphs(G))[0]
        #print (len(g.nodes()), 'nodes')
        #print (len(g.edges()), 'edges')
        return(g)
    else:
        #print (len(G.nodes()), 'nodes')
        #print (len(G.edges()), 'edges')
        return(G)


def calculate_distances (G, nodes_from, nodes_to, 
                         sp=None, node2index=None, lcc=True):
    
    """
    pair of nodes that do not have a path
    do not contribute to the final value
    
    sp: numpy matrix
    index2node: dict
    
    """
    
    ds = defaultdict(dict)
    

    for i in nodes_from:
        for j in nodes_to:
            if i == j:
                ds[i][j] = 0
            else:
                if nx.has_path(G,i, j):
                    if sp is None:
                        ds[i][j] = nx.shortest_path_length(G,i, j)
                    else:
                        ds[i][j] = sp[node2index[i],node2index[j]]
                else:
                    ds[i][j] = float('nan')
        
    ds = pd.DataFrame.from_dict(ds)
    # nodes_to: rows
    # nodes_from: columns 
    
    dic = {}
    
    dic['shortest'] = ds.mean().mean()
    dic['closest'] = ds.min().mean()
    
    return (dic)

def calculate_significance_pars(network, nodes, nodes_random=None, bins=None, n_random=1000, 
                               min_bin_size=100, seed=452456):
    
    
    if nodes_random is None:
        network_nodes = list(network.nodes())
        nodes_random = []
        for i in range(n_random):
            shuffle(network_nodes)
            nodes_random.append(network_nodes[:len(nodes)])
    res = {}
    ## lcc size
    network_sub = network.subgraph(nodes)
    component_nodes = network_utilities.get_connected_components(network_sub, False)[0]
    res['lcc'] = len(component_nodes)
    
    ## proportion nodes in lcc
    res['lcc_p'] = len(set(nodes) & set(component_nodes))/len(component_nodes)
    
    ## conductance
    res['cond'] = nx.algorithms.cuts.conductance(network, nodes)
    
    ## density
    res['density'] = nx.density(nx.subgraph(network, nodes))
    
    values = {}
    ## vectors for output
    values['lcc'] = np.empty(len(nodes_random))
    values['lcc_p'] = np.empty(len(nodes_random))
    values['cond'] = np.empty(len(nodes_random))
    values['density'] = np.empty(len(nodes_random))
    
    ## loop in random
    for i, nodes in enumerate(nodes_random):
        network_sub = network.subgraph(nodes)
        component_nodes = network_utilities.get_connected_components(network_sub, False)[0]
        values['lcc'][i] = len(component_nodes)
        values['lcc_p'][i] = len(set(nodes) & set(component_nodes))/len(component_nodes)
        values['cond'][i] = nx.algorithms.cuts.conductance(network, nodes)
        values['density'][i] = nx.density(network_sub)
    
    final = defaultdict(dict)
    for par in res.keys():
        m, s = np.mean(values[par]), np.std(values[par])
        if s == 0:
            z = 0.0
        else:
            z = (res[par] - m) / s
        final[par]['mean'], final[par]['std'], final[par]['z'] = m,s,z
        final[par]['value'] = res[par]
    
    return (final) 


def calculate_sab(G, nodes_from, nodes_to):

    # distances WITHIN the two gene sets:
    d_A = separation.calc_single_set_distance(G,set(nodes_from))
    d_B = separation.calc_single_set_distance(G,set(nodes_to))

    # distances BETWEEN the two gene sets:
    d_AB = separation.calc_set_pair_distances(G,set(nodes_from),set(nodes_to))

    # calculate separation
    s_AB = d_AB - (d_A + d_B)/2.
    
    return(s_AB)


def calculate_sab_null(G, nodes_from, nodes_to,n_random=1000):
    
    all_genes = list(G.nodes())
    null_ref = []
    for i in range(n_random):
        seeds = set(random.sample(all_genes,len(nodes_from)))
        targets = set(random.sample(all_genes, len(nodes_to)))
        sab = calculate_sab(G, seeds, targets)
        null_ref.append(sab)
    
    return (null_ref)


def ref_distribution(S,T,G,cpu=2, iterations=10):
    
    """
    S: [list] set of disease proteins
    T: [list] set of drug targets
    G: [nx.Graph] interactome
    reference distribution to assess the significance between (S,T)
    """
    
    iterable = [(S,G) for i in range(iterations)]
    p = Pool(cpu)
    nS_v = p.map(get_random_nodes, iterable)
    #logging.info('Random disease genes - %d cpu %d iterations'%(cpu,iterations))
    #nS_v = list(tqdm.tqdm(p.imap(get_random_nodes, iterable), total=len(iterable)))
    p.close()
    
    #logging.info('Random target genes - %d cpu %d iterations'%(cpu,iterations))
    iterable = [(T,G) for i in range(iterations)]
    p = Pool(cpu)
    nT_v = p.map(get_random_nodes, iterable)
    #nT_v = list(tqdm.tqdm(p.imap(get_random_nodes, iterable), total=len(iterable)))
    p.close()
    
    ref_s = []
    ref_c = []
    ref_k = []
    for i in range(iterations):
        nS = nS_v[i]
        nT = nT_v[i]
        ds, dc, dk = compute_distances(nS,nT, G)
        ref_s.append(ds)
        ref_c.append(dc)
        ref_k.append(dk)
    
    
    return([ref_s, ref_c, ref_k])
    
    
def compute_sp_to_closest (S,G):
    """
    S: [list] set of source nodes
    G: [nx.Graph] interactome
    """
    if len(S) == 1:
        raise NetworkUtilsInputError('Input received just one node')
    
    
    dt = pd.DataFrame(index=S, columns=S)
    
    
    for s in dt.index:
        for t in dt.columns:
            d = nx.shortest_path_length(G,source=s, target=t)
            dt[t].loc[s] = d
    
    ## ignoring diagonal
    m = np.asmatrix(dt)
    np.fill_diagonal(m, np.inf)
    dt = pd.DataFrame(m,index=dt.index, columns=dt.columns)
    
    
    dc = dt.apply(min)
    
    return (dc)



def get_lcc(G,S):
    """
    S: [list] set of source nodes
    G: [nx.Graph] interactome
    """
    if len(S) == 0:
        return (nx.Graph())
    else:
        g = nx.subgraph(G,S)
        if len(g.nodes()) > 0:
            lcc = max(nx.connected_component_subgraphs(g), key=len)
            return (lcc)
        else:
            return(g)

        
        
def get_lcc_significance(G,seeds,n_random=1000):

    # getting all genes in the network  
    all_genes = G.nodes()

    number_of_seed_genes = len(seeds & set(all_genes))
    
    l_list  = []

    # simulations with randomly distributed seed nodes
    for i in range(1,n_random+1):

        # get random seeds
        rand_seeds = set(random.sample(all_genes,number_of_seed_genes))

        # get rand lcc
        lcc = get_lcc(G,rand_seeds)
        lcc = len(lcc.nodes())
        l_list.append(lcc) 
        

    # get the actual value
    lcc_observed = get_lcc(G,seeds)
    lcc_observed_size = len(lcc_observed.nodes())

    # get the lcc z-score:
    l_mean = np.mean(l_list)
    l_std  = np.std(l_list)

    if l_std == 0:
        z_score = float('nan')
    else:
        z_score = (1.*lcc_observed_size - l_mean)/l_std

    
    return ({'lcc_size':lcc_observed_size, 'z_score':z_score})


    
def read_edgelist (infile,sep=' ',header=False):
    if header:
        lines = open(infile, 'r').readlines()[1:]
    else:
        lines = open(infile, 'r').readlines()
    edges = []
    for line in lines:
        a, b = line.rstrip().split(sep)
        edges.append((a,b))
    g = nx.from_edgelist(edges)
    return(g)


def get_visualization_subnetwork(seeds,targets,G):
#if True:
    #seeds = list(set(disease2genes['gastroenteritis']) & set(G.nodes()))
    #targets = list(set(chemical2genes['genistein']) & set(G.nodes()))
    ### caculate distances
    dist = defaultdict(dict)
    for i in targets:
        for j in seeds:
            dist[i][j] = nx.shortest_path_length(G,i,j)

    dist = pd.DataFrame.from_dict(dist)
    
    ### for each target, get the shortest path to the closest disease protein
    edges = []
    for i in dist.columns:
        for j in dist.index:
            if dist[i].loc[j] == dist.loc[j].min():
                p = nx.shortest_path(G, i, j)
                for k in range(len(p)-1):
                    edges.append((p[k],p[k+1]))
                    
    ## include edges in LCCs
    lccd = get_lcc(G, seeds)
    edges = edges + list(lccd.edges())
    
    lcct = get_lcc(G, targets)
    edges = edges + list(lcct.edges())
    
    
    edges = list(set(edges))
    g = nx.Graph()
    g.add_edges_from(edges)
    
    ## create edgelist table
    dic = defaultdict(dict)
    c = 0
    for edge in edges:
        dic[c]['source'] = int(edge[0])
        dic[c]['target'] = int(edge[1])
        dic[c]['intra_disease'] = edge[0] in seeds and edge[1] in seeds
        dic[c]['intra_targets'] = edge[0] in targets and edge[1] in targets
        dic[c]['inter'] = 0
        if edge[0] in seeds and edge[1] in targets:
            dic[c]['inter'] = 1
        if edge[1] in seeds and edge[0] in targets:
            dic[c]['inter'] = 1
        c = c + 1
    edgetable = pd.DataFrame.from_dict(dic, orient='index')
    
    
    
    ## create node table
    ntable = defaultdict(dict)
    c = 0
    mapping = convert_gene_id(list(g.nodes()), 
                                              in_format='entrez', 
                                              out_format='symbol')
    for i in g.nodes():
        ntable[c]['node'] = int(i)
        ntable[c]['degree_sub'] = g.degree(i)
        ntable[c]['degree_orig'] = G.degree(i)
        ntable[c]['disease'] = i in seeds
        ntable[c]['target'] = i in targets
        ntable[c]['symbol'] = float('nan')
        if str(int(i)) in mapping.keys():
            ntable[c]['symbol'] = mapping[str(int(i))]
        c = c + 1
    
    
    nodetable = pd.DataFrame.from_dict(ntable, orient='index')
    
    return(edgetable, nodetable)

def convert_gene_api(query, out_format, mg):

    out = None
    if out_format == 'entrez':
        out_format = 'entrezgene'
    res = mg.query(query)
    if len(res['hits']) > 0:
        for h in res['hits']:
            if h['taxid'] == 9606 and out_format in h.keys():
                out = h[out_format]
                
    return(out)
        
def convert_gene_id(query, in_format='symbol', 
                    out_format='entrez'):
    
    """
    in_format: [symbol, entrez, ensemblg, ensemblt, ensemblp]
    """
    
    dt = pd.read_csv(infolder + '/databases/geneids.csv')
    
    dt = dt[[in_format, out_format]]
    dt = dt[~dt.isnull().any(axis=1)]
    
    dt['entrez'] = [str(int(i)) for i in dt['entrez']]
    mapping = {i:j for i,j in zip(dt[in_format], dt[out_format])}
    
    
    print ('Internal database')
    pbar = ProgressBar()
    res = {}
    missing = []
    if type(query) == list:
        query = list(map(str, query))
        for gene in pbar(query):
            if gene in mapping.keys():
                res[gene] = mapping[gene]
            else:
                missing.append(gene)
    else:
        query = str(query)
        if query in pbar(mapping.keys()):
            res[query] = mapping[query]
        else:
            missing.append(query)
    
    
    print ('Found %d in internal database'%len(res))
    
    print ('API')
    c = 0
    if len(missing) > 0:
        pbar = ProgressBar()
        mg = mygene.MyGeneInfo()
        for elem in pbar(missing):
            f = convert_gene_api(elem, out_format, mg)
            if f:
                res[elem] = f
                c = c + 1
    
    print ('Found %d in API'%c)
    
    return(res)