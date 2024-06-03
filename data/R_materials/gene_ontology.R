library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)

# only run this code if you get a HDO.db related error while running library(clusterProfiler)
BiocManager::install("HDO.db")

# Arabidopsis Database
AT_DB = "org.At.tair.db"

# use argument force=TRUE if warning message arrives
BiocManager::install(AT_DB, character.only = TRUE)
library(AT_DB, character.only = TRUE)
keytypes(org.At.tair.db)

# reading in data from deseq2
setwd("~/Downloads/05_counts/")
res = read.csv("star_featurecounts_results.csv", header=TRUE, row.names = 1)

# rename the gene names (row names here) so they match the TAIR database
rownames(res) <- substr(rownames(res), 6, nchar(rownames(res)))
head(res)

# we want the log2 fold change 
original_gene_list <- res$log2FoldChange

# name the vector
names(original_gene_list) <- row.names(res)
original_gene_list

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in descending order
gene_list = sort(gene_list, decreasing = TRUE)
df <- as.data.frame(names(gene_list))
df
write.table(df, file="gene_list.csv", row.names = FALSE, col.names = FALSE)

gsea <- gseGO(geneList=gene_list, 
                       ont ="ALL", 
                       keyType = "TAIR",  
                       minGSSize = 3, 
                       maxGSSize = 1000, 
                       pvalueCutoff = 0.05, 
                       verbose = TRUE, 
                       OrgDb = org.At.tair.db)
summary(gsea)
dotplot(gsea, showCategory=10, split=".sign") + facet_grid(.~.sign)

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "ath",
                     minGSSize = 3,
                     maxGSSize = 1000,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)
summary(gsea_kegg)
dotplot(gsea_kegg, showCategory = 10, split=".sign") + facet_grid(.~.sign)

# now let's do an over-representation analysis (gene names only)
# grab the significant DGEs only
res_sig = subset(res, abs(res$log2FoldChange) > 1)
res_sig = subset(res_sig, res_sig$padj < 0.05)
gene <- rownames(res_sig)

# handy link to organisms in KEGG db https://www.genome.jp/kegg/catalog/org_list.html
ora <- enrichKEGG(gene,
                 organism="ath",
                 pvalueCutoff=0.05,
                 pAdjustMethod="fdr")

browseKEGG(ora, 'ath04075')
browseKEGG(ora, 'ath04626')
browseKEGG(ora, 'ath00906')
browseKEGG(ora, 'ath04016')
browseKEGG(ora, 'ath04712')
browseKEGG(ora, 'ath00270')

