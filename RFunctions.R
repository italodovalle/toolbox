library(clusterProfiler)
library("org.Hs.eg.db)
library(ReactomePA)

data(geneList, package="DOSE")

enrichment_test_GO <- function(genelist, adjustmethod='BH',pcutoff=0.05,qcutoff=0.05,
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


          ego <- enrichGO(gene          = genes,
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


run_enrichment <- function(drug){


        outdir = '../output/broad_screening/'
        ont = 'BP'
        entrez_ids = as.character(dx[i, 'target_entrez'])
        entrez_ids = unlist(strsplit(entrez_ids, '\\|'))
        ego = enrichment_test_GO(entrez_ids,
                             "BH",0.05,0.05,
                             ont = ont,input="ENTREZ",
                             readable = TRUE)
        outname = paste(drug, ont, sep = '_')
        outname = paste0(outname, '.csv')
        outname = paste0(outdir, outname)
        write.csv(ego, outname)

}

enrichment_test_Reactome <- function(gene,pcutoff=0.05,adjustmethod="BH",qcutoff=0.05,
                                       min_gs_size=10, max_gs_size=500, input = 'ENTREZ') {

        if (input == "SYMBOL"){
            gene.df <- bitr(genelist, fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)
            genes = gene.df$ENTREZID
         }
         else {
            genes = genelist

         }


        x <- enrichPathway(genes,pvalueCutoff=pcutoff,
                            pAdjustMethod = adjustmethod,
                            qvalueCutoff = qcutoff,
                            minGSSize = min_gs_size,
                            maxGSSize = max_gs_size,
                            organism = "human",
                            readable=T)
        df = as.data.frame(x)
        return(df)

}
