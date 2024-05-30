library(DESeq2)
library(pheatmap)


setwd("~/Downloads/05_counts/")

# open the featureCounts matrix as a matrix
counts <- read.table("combined.counts.txt", sep="\t", header=TRUE, row.names="Geneid")

# keep only the sample count columns
counts <- counts[, (ncol(counts)-5):ncol(counts)]

sample_id <- colnames(counts)
sample_id <- sapply(strsplit(sample_id, "[.]"), "[", 1)
colnames(counts) <- sample_id
counts <- as.matrix(counts)

# create a sample dataframe where sample ids are in the SAME ORDER as counts matrix above
hours_cold <- sapply(strsplit(sample_id, "_"), "[", 2)
rep <- sapply(strsplit(sample_id, "_"), "[", 3)
samples <- data.frame(sample_id, hours_cold, rep)
rownames(samples) <- samples$sample_id
samples <- samples[, -1]

# create a DESeq object from our matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ hours_cold)

# run the DESeq!
dds <- DESeq(dds)
plotDispEsts(dds)

rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = c("hours_cold"))

res <- results(dds, alpha = 0.05, contrast = c("hours_cold", "0h", "3h"))
res_ord <- res[order(res$padj, -abs(res$log2FoldChange)),]
write.csv(res_ord, file="star_featurecounts_results.csv")
summary(res)
plotMA(res, ylim=c(-10, 10))

res_sig <- as.data.frame(res[ which(res$padj < 0.05),])
res_sig <- res_sig[order(res_sig$padj, -abs(res_sig$log2FoldChange)),]

write.csv(res_sig, file="star_featurecounts_sig_results.csv")
summary(res_sig)

# these steps can be used to find individual points on the graph
# after running "identify", click on the plot, then hit "finish" button in top right of plot
plotMA(res, ylim=c(-10, 10))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


# optional volcano plot
library(apeglm)
library(EnhancedVolcano)
library(ggpubr)

EnhancedVolcano(res, x = 'log2FoldChange',
                lab = row.names(res),
                pCutoff = 1e-100,
                FCcutoff = 2,
                y = 'padj',)

# create a plot for a single gene
plotCounts(dds, gene="gene-AT4G25480", intgroup="hours_cold")

# create a heatmap
# first run vst (quickly estimate dispersion trend and apply a variance stabilizing transformation)
vsd <- vst(dds)

mat <- assay(vsd)[ head(order(-abs(res$log2FoldChange)), 20), ]
mat <- mat - rowMeans(mat) 
mat <- data.frame(mat)

pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)
