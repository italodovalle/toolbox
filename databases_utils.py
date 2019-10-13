#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 17:31:43 2017

@author: italodovalle
"""

import os
import json
import mygene
import itertools
import pandas as pd
import pubchempy as pcp
import networkx as nx
from collections import defaultdict
from progressbar import ProgressBar
from mygene import MyGeneInfo


infolder = os.getcwd() + '/toolbox'


def convert_gene_api(query):
    
    mg = MyGeneInfo()
    dic = {}
    out = float('nan')
    out_format = 'entrezgene'
    try:
        res = mg.query(query)
    except:
        res = {}
        res['hits'] = []
    if len(res['hits']) > 0:
        for h in res['hits']:
            if h['taxid'] == 9606 and out_format in h.keys():
                out = h[out_format]
    else:
        out = float('nan')
        
    dic[query] = out
                
    return (dic)

def uniprot2entrez(uniprot_list, mapping_file = None):
    
    
    if mapping_file is None:
        df = pd.read_csv(infolder + 'entrez_uniprot_mapping.csv', 
                         index_col = 0)
    else:
        df = pd.read_csv(mapping_file, index_col = 0)
    
    dx = pd.DataFrame()
    dx['uniprot'] = uniprot_list
    dx['original'] = 1
    res = pd.merge(dx, df, on= 'uniprot', how='outer')
    res = res[res['original'] == 1]
    res = res[['uniprot', 'entrez']].drop_duplicates()
    mapping = {i:j for i,j in zip(res['uniprot'], res['entrez'])}
    
    #print (len(mapping))
    
    c = 0
    print ('Complete with API')
    if res[res.entrez.isnull()].shape[0] > 0:
        pbar = ProgressBar()
        for g in pbar(res[res.entrez.isnull()]['uniprot']):
            x = convert_gene_api(g)
            mapping.update(x)
            if x[g]:
                c = c + 1
    
    
    return (mapping)








def query_pubchem (query):
    cid = float('nan')
    try:
        c = pcp.get_compounds(query, 'name')
    except:
        c = []
    if len(c)>0:
        cid = c[0].cid
    return(cid)


def load_mesh_graph(infile):
    disease2code = defaultdict(list)
    code2name = {}
    codes = []
    for i in open(infile).readlines():
        name, code = i.rstrip().split(';')
        if code.startswith('C'):
            disease2code[name.lower()].append(code)
            code2name[code] = name.lower()
            codes.append(code)

                 
    edges = []
    for i in set(codes):
        if len(i) > 4:
            a, b = i[:-4], i
            edges.append((a,b))
            
    g = nx.DiGraph()
    g.add_edges_from(edges)  
    
    return(g, disease2code, code2name)
    

### old
def set_db_folder(path):
    """
    it must follow a particular folder structure
    """
    global __DBPATH__
    if path.endswith('/'):
        __DBPATH__ = path
    else:
        __DBPATH__ = path + '/'
        
        
def load_ncbi():
    with open (__DBPATH__ + 'ncbi/ncbi2symbol.json') as fp:
        global ncbi2symbol
        ncbi2symbol = json.load(fp)
    with open(__DBPATH__ + 'ncbi/old2new.json') as fp:
        global old2new
        old2new = json.load(fp)
        

        
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
        
    
    
    


def get_symbol(entrez):
    entrez = str(entrez)
    if entrez in ncbi2symbol.keys():
        return(ncbi2symbol[entrez])
    elif entrez in old2new.keys():
        if old2new[entrez] != '-':
            newentrez = old2new[entrez]
            return(ncbi2symbol[newentrez])   
            


def get_entrezid(symbols):
    
    mg = mygene.MyGeneInfo()
    if type(symbols) == str:
        symbols = [symbols]
    
    entrez = {}
    
    for symbol in symbols:
        try:
            entrez[symbol] = int(mg.query(symbol)['hits'][0]['_id'])
        except:
            entrez[symbol] = None
            
    return (entrez)
    


def kegg_compound2genes(compound):
    
    with open(__DBPATH__ + 'kegg/enzymes_compounds.json', 'r') as fp:
        kegg = json.load(fp)
    target = []
    for i in kegg.keys():
        if compound in kegg[i]:
            target.append(i)
    target = list(set(target))
    
    return(target)    
    
    
    
def get_go_similarity(S,T):
    
    def s_go(a,b):
        if a not in gene_go.keys() or b not in gene_go.keys():
            shared = []
        else:
            shared = list(set(gene_go[a]) & set(gene_go[b]))
        if len(shared) == 0:
            return (0)
        else:
            sizes = [len(go_gene[x]) for x in shared]
            similarity = 2./min(sizes)
            return (similarity)
    
    
    with open(__DBPATH__ + 'go/gene2go.json','r') as fp:
        gene_go = json.load(fp)
    
    with open(__DBPATH__ + 'go/go2gene.json','r') as fp:
        go_gene = json.load(fp) 


    pairs = list(itertools.product(T,S))
    s = []
    for pair in pairs:
        if pair[0] != pair[1]:
            s.append(s_go(pair[0],pair[1]))
    #s_S_T = sum(s)/len(pairs)
    
    return(s)
    

def get_recon_species_dic(infile):
    recon = pd.read_csv(infile,index_col = 0)
    ids = [i[:-2] for i in recon.index]
    recon2name = dict(zip(ids, recon['name']))
    return(recon2name)
    
