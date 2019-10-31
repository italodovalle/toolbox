#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import random
from random import shuffle
#import separation
from guney_code.wrappers import network_utilities, get_random_nodes, calculate_closest_distance, calculate_separation_proximity, calculate_proximity
from scipy import stats
from scipy import sparse
import mygene
import os
from progressbar import ProgressBar
import pickle


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
