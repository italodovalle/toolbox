#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 14:11:00 2018

@author: italodovalle
"""

import rpy2
import os
#os.environ['R_HOME'] = '/usr/local/Cellar/r/3.4.1_1/lib/R'
#os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
from rpy2 import robjects
from rpy2.robjects import pandas2ri


robjects.r('''
        library(ReactomePA)
        library(clusterProfiler)
        library("org.Hs.eg.db")

        data(geneList, package="DOSE")

        enrichment_test_Reactome <- function(gene,pcutoff,adjustmethod,qcutoff,
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


    enrichment_test_GO <- function(genelist, adjustmethod,pcutoff,qcutoff,
                   ont = "BP",input="SYMBOL",
                   readable = TRUE){


         if (input == "SYMBOL"){
            gene.df <- bitr(genelist, fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)
            genes = gene.df$ENTREZID
         }
         if (input == "ENTREZ"){
            genes = genelist

         }


          ego <- enrichGO(gene          = gene.df$ENTREZID,
                          universe      = names(geneList),
                          OrgDb         = org.Hs.eg.db,
                          ont           = ont,
                          pAdjustMethod = adjustmethod,
                          pvalueCutoff  = pcutoff,
                          qvalueCutoff  = qcutoff,
                          readable      = readable)


         ego = as.data.frame(ego)
         return (ego)

}
''')


def get_enrichment_Reactome (input_vector,pcutoff=0.05,adjustmethod="BH",
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

    enrich = robjects.r['enrichment_test_Reactome']
    genes = robjects.StrVector(input_vector)
    enrichment = pandas2ri.ri2py_dataframe(enrich(genes,
                                                  pcutoff,adjustmethod,
                                                  qcutoff,
                                                  min_gs_size = 10,
                                                  max_gs_size=500))
    return(enrichment)



def get_enrichment_GO (input_vector,pcutoff=0.05,adjustmethod="BH",
                    qcutoff=0.2,ont = "BP",input="SYMBOL",
                    readable = True):

    """
    Gene Ontology Enrichment Analysis (clusterProfiler Bioconductor)

    Args:
        input_vector: (obj:list) gene ids str format
        pcutoff: p-value threshold
        adjustmethod: Multiple Testing Correction method
                      one of "holm", "hochberg", "hommel", "bonferroni",
                      "BH", "BY", "fdr", "none"
        input: 'SYMBOL' or 'ENTREZ'
        ont: Gene Ontology Category
             "BP": biological process
             "CC": cellular compartment
             "MF: molecular function

    Returns:
        df: DataFrame with enrichment results
    """

    enrich = robjects.r['enrichment_test_GO']
    genes = robjects.StrVector(input_vector)
    enrichment = pandas2ri.ri2py_dataframe(enrich(genes,
                                                  adjustmethod,
                                                  pcutoff,
                                                  qcutoff,
                                                  ont=ont,
                                                  input=input,
                                                  readable=readable))
    return(enrichment)
