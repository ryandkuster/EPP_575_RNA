library(tidyverse)
library(tximport)
library(GenomicFeatures)
library(pheatmap)
library(DESeq2)

# create a workding directory
setwd("/Users/ryankuster/Documents/projects/teaching/2024_06_EPP_575_RNA/data/derived")

# give the full path to the "salmon_results" (it should end in "salmon_results"
dir <- "/Users/ryankuster/Documents/projects/teaching/2024_06_EPP_575_RNA/data/derived/salmon_results"

# load the gff3 file, then create a transcript database/dataframe for use with deseq
txdb <- makeTxDbFromGFF("/Users/ryankuster/Documents/projects/teaching/2024_06_EPP_575_RNA/data/reference/genomic.gff")
keytypes(txdb)
k <- keys(txdb, keytype = "CDSNAME")
str(k)

txdf = AnnotationDbi::select(txdb, k, "GENEID", "CDSNAME")

# load in the metadata
samples <- read_csv("/Users/ryankuster/Documents/projects/teaching/2024_06_EPP_575_RNA/data/derived/salmon_results/rename_samples_file.csv")

Qfiles <- file.path(dir, samples$quant_file)

# this step imports the count data from salmon
txi <- tximport(files = Qfiles, type = "salmon", tx2gene = txdf)
colnames(txi$counts) <- samples$sample_id
names(txi)
head(txi$counts)
summary(txi)

# now we convert the txi object into a deseq-formatted object
dds_data <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~condition)
dds <- DESeq(dds_data)

# plot dispersion
plotDispEsts(dds)

# summarize results
res <- results(dds)
head(res)

# create a contrast with lfcThreshold and alpha cutoff (first list item is condition from samples object)
res_sig <- results(dds, alpha = 0.05, contrast = c("condition", "Col_0h", "Col_3h"))

max_res_sig <- as.data.frame(res_sig[ which(res_sig$padj < 0.05 & abs(res_sig$log2FoldChange) > 2),])
write.table(max_res_sig, file=paste0(dir,"/salmon_sig_results.csv"), row.names = T)
summary(res_sig)
plotMA(res_sig, ylim=c(-12,12))


