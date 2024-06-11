---
title: "GO_enrichment_genes_with_SVs"
author: "Beant Kapoor"
date: "2023-06-29"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r Install and load required packages}
# Install
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("AnnotationForge")
#BiocManager::install("GO.db")

#BiocManager::install("clusterProfiler")
#install('RSQLite', force = TRUE)

#BiocManager::install("enrichplot", force = TRUE)

# load
library(AnnotationForge)
library(GO.db)
library(dplyr)
library(splitstackshape)
library(clusterProfiler)
library(RSQLite)
library(enrichplot)
library(ggplot2)
```

```{r Quercus rubra gene annotation custom database}
##########################################
# IF RUN ONCE THEN NO NEED TO RUN THIS CHUNK AGAIN
##########################################

# read in file
annot <- read.table("Qrubra_687_v2.1.annotation_info.txt", sep = "\t", header = FALSE, fill = TRUE)

# extract V1, V2, V12 and put them in separate df
rubra_genes <- annot[, c(1,2,12)]

# remove if blanks are present in 3rd column
rubra_genes <- rubra_genes[(rubra_genes$V12 != ""), ]

# remove duplicated genes
rubra_genes <- rubra_genes[order(rubra_genes$V2), ]
rubra_genes <- rubra_genes[!duplicated(rubra_genes$V2), ]

# provide right colnames, this is very important
colnames(rubra_genes) <- c("GID","SYMBOL","GENENAME")

# reset row names
rownames(rubra_genes) <- NULL #rubra_genes is ready

# extract V1 and V2 columns from annot
rubra_genes_chr <- annot[, c(1,2)]

# sub string the 2 column to get chr number only
rubra_genes_chr$V2 <- substr(rubra_genes_chr$V2,7,8)

# remove duplicated entries in 1 column
rubra_genes_chr <- rubra_genes_chr[order(rubra_genes_chr$V1), ]
rubra_genes_chr <- rubra_genes_chr[!duplicated(rubra_genes_chr$V1), ]

# just keep those genes which are in chromosomes and not unplaced scaffolds
colnames(rubra_genes_chr) <- c("GID","CHROMOSOME")
rubra_genes_chr <- rubra_genes_chr %>% filter(grepl('01|02|03|04|05|06|07|08|09|10|11|12', CHROMOSOME))

# rubra_genes_chr is ready

# extract V1 and V10
rubra_go <- annot[,c(1,10)]

# remove rows with empty cells in V2 column
rubra_go <- rubra_go[(rubra_go$V10 != ""), ]
colnames(rubra_go) <- c("GID","GO")

# split GO terms to different rows
rubra_go <- cSplit(rubra_go, "GO", " ", "long")

# remove rows with NAs
rubra_go <- rubra_go[!is.na(rubra_go$GID), ]

# add another column EVIDENCE
rubra_go <- rubra_go %>% mutate(EVIDENCE = "IEA")

# check entire column of GO IDS
any(!grepl("^GO:", as.character(rubra_go$GO)))

# filter out any GO IDs which does not start with GO
rubra_go <- rubra_go[grepl("^GO:", as.character(rubra_go$GO)),]

# rubra_go is ready

# generate custom database
makeOrgPackage(gene_info=rubra_genes, chromosome=rubra_genes_chr, go=rubra_go,
               version="2.1",
               maintainer="Beant Kapoor <bkapoor@vols.utk.edu>",
               author="Beant Kapoor <bkapoor@vols.utk.edu>",
               outputDir = ".",
               tax_id="3512",
               genus="Quercus",
               species="rubra",
               goTable="go",
               verbose=TRUE)
```

```{r Genes with inversion breakpoints}
install.packages("./org.Qrubra.eg.db", repo = NULL, type = "source")

library(org.Qrubra.eg.db)

genes_having_inversions <- read.table("gene_names_having_inversion_breakpoints.txt", header = FALSE)

GO_1 <- enrichGO(gene = genes_having_inversions$V1,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'ALL',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05,
                 readable = TRUE)

#goplot(GO_1)

tiff("GO term enrichment genes having inversion breakpoints.tiff", units = "in", width = 9, height = 7, res = 300)
dotplot(GO_1) + ggtitle("Genes having inversion breakpoints")
dev.off()
```

```{r Genes with insertions}
genes_having_insertions <- read.table("gene_names_having_insertions.txt", header = FALSE)

GO_2 <- enrichGO(gene = genes_having_insertions$V1,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'ALL',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05)

tiff("GO term enrichment genes having insertions.tiff", units = "in", width = 9, height = 7, res = 300)
dotplot(GO_2) + ggtitle("Genes having insertions")
dev.off()
```

```{r Genes with breakend breakpoints}
genes_having_breakend <- read.table("gene_names_having_breakend_breakpoints.txt", header = FALSE)

GO_3 <- enrichGO(gene = genes_having_breakend$V1,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'ALL',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05)

dotplot(GO_3) + ggtitle("Genes having breakend mates")

#####################################
# NO GO TERMS ARE ENRICHED IN THE GENES AFFECTED BY BREAKENDS
#####################################
```

```{r Genes with ALL SVs}
genes_having_ALL_SVs <- read.table("gene_names_having_ALL_SVs.txt", header = FALSE)

GO_4 <- enrichGO(gene = genes_having_ALL_SVs$V1,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'BP',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'fdr')

# number of enriched GO terms
length(GO_4$ID)

tiff("GO term enrichment genes having ALL SVs.tiff", units = "in", width = 9, height = 7, res = 300)
dotplot(GO_4) + ggtitle("All genes affected by SVs")
dev.off()

png('GO term enrichment genes having ALL SVs.png', width = 10, height = 10, units = "in", res = 300)
dotplot(GO_4) + ggtitle("All genes affected by SVs")
dev.off()
```

```{r Genes not affected by SVs}
#genes_NOT_affected_SVs <- as.data.frame(setdiff(rubra_genes$V2, genes_having_ALL_SVs$V1))
#names(genes_NOT_affected_SVs) <- "gene"

# save this dataframe so that we don't have to run the above code again
#write.table(genes_NOT_affected_SVs, "gene_names_NO_SVs.txt", row.names = FALSE, quote = FALSE)

# read in data
genes_having_NO_SVs <- read.table("gene_names_NO_SVs.txt", header = TRUE)

GO_5 <- enrichGO(gene = genes_having_NO_SVs$gene,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'BP',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'fdr')

# tiff
tiff("GO term enrichment genes NO SVs.tiff", units = "in", width = 9, height = 7, res = 300)
dotplot(GO_5) + ggtitle("Genes not affected by SVs")
dev.off()

# png
png('GO term enrichment genes NO SVs.png', width = 10, height = 10, units = "in", res = 300)
dotplot(GO_5) + ggtitle("Genes not affected by SVs")
dev.off()

length(GO_5$ID)
```

```{r genes inside SV hotspots}
genes_inside_SV_hotspots <- read.table("genes_inside_SV_hotspots.txt", header = FALSE)

GO_6 <- enrichGO(gene = genes_inside_SV_hotspots$V1,
                 OrgDb = org.Qrubra.eg.db,
                 keyType = 'SYMBOL',
                 minGSSize = 1,
                 ont = 'BP',
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'fdr')

tiff("GO term enrichment genes inside SV hotspots.tiff", units = "in", width = 9, height = 7, res = 300)
dotplot(GO_6) + ggtitle("Genes inside SV hotspots")
dev.off()
```
