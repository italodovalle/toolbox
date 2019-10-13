#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:11:00 2018

@author: italodovalle
"""

import rpy2
import os
os.environ['R_HOME'] = '/usr/local/Cellar/r/3.4.1_1/lib/R'
from rpy2 import robjects
from rpy2.robjects import pandas2ri


robjects.r('''               
        library(ReactomePA)
        enrichment_test <- function(gene,pcutoff,adjustmethod,qcutoff,
                       min_gs_size, max_gs_size) {
                      x <- enrichPathway(gene,pvalueCutoff=pcutoff,
                                         pAdjustMethod = adjustmethod,
                                         qvalueCutoff = qcutoff,
                                         minGSSize = min_gs_size,
                                         maxGSSize = max_gs_size,
                                         organism = "human",
                                         readable=T)
                      df = as.data.frame(x)
                      return(df)
                      
        }
        ''')


def get_enrichment (input_vector,pcutoff=0.05,adjustmethod="BH",
                    qcutoff=0.2, min_gs_size = 10, max_gs_size=500,
                    organism="human"):
    
    """
    Reactome Enrichment Analysis (ReactomePA Bioconductor)
    
    Args:
        input_vector: (obj:list) gene entrez ids str format
        pcutoff: p-value threshold
        adjustmethod: Multiple Testing Correction method
                      one of "holm", "hochberg", "hommel", "bonferroni", 
                      "BH", "BY", "fdr", "none"
                      
    Returns:
        df: DataFrame with enrichment results
    """
    
    enrich = robjects.r['enrichment_test']
    genes = robjects.StrVector(input_vector)
    enrichment = pandas2ri.ri2py_dataframe(enrich(genes,
                                                  pcutoff,adjustmethod,
                                                  qcutoff,
                                                  min_gs_size = 10, 
                                                  max_gs_size=500))
    return(enrichment)
