rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)
library(stringr)
library(ramify) # this is for clip, for the volcano

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/akthelia/results/000_quantification"
results_dir = '/Users/adrian/research/akthelia/results/002_DEGs'

#
# 1. generate gene to transcript mapping
#
t2g = read.csv('/Users/adrian/software/kallisto/human_index_standard/annotation.tsv', sep='\t')
a = do.call(rbind, stringr::str_split(t2g$description, '\\[Source'))
t2g$description2 = a[, 1]
dim(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
paths = file.path(dirnames, 'kallisto_output_a/abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[8]))
labels = str_remove(labels, '_processed')

replicates = rep(c('A', 'B', 'C'), 8)
genotypes = c(rep('T0570', 12), rep('T84', 12))
times = rep(c(rep('24', 3), rep('06', 3)), 4)
treatments = rep(c(rep('011', 6), rep('DMSO', 6)), 2)
metadata = data.frame(labels)
metadata$path = paths
metadata$replicate = replicates
metadata$genotype = genotypes
metadata$time = times
metadata$treatment = treatments
dim(metadata)

# get the appropriate samples
metadata = metadata[metadata$genotype == 'T0570', ]
metadata = metadata[metadata$time == '24', ]
dim(metadata)
View(metadata)

clustera_indexes = 1:3
clusterb_indexes = 4:6

#
# 3. contrasts
#
count_threshold = 20
effect_size_threshold = log2(2)
tpm_threshold = 2

#
# 3.1. contrast 
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g[, 2:3], ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=metadata, design=~treatment) 
dds$treatment = relevel(dds$treatment, ref="DMSO")

# keep features with at least 20 counts median difference
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
a = counts(dds)[ , clustera_indexes]
b = counts(dds)[ , clusterb_indexes]
c = rowMedians(a) - rowMedians(b)
keep = abs(c) >= count_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

# keep features with at least a max median expression of 1 TPM.
cat(blue(paste('size before counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)
subset = txi$abundance[names(dds), ]
a = rowMedians(subset[ , clustera_indexes])
b = rowMedians(subset[ , clusterb_indexes])
c = pmax(a, b)
keep = c >= tpm_threshold
dds = dds[keep, ]
cat(blue(paste('size after counts filtering:', dim(dds)[1], sep=' ')), fill=TRUE)

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) # it does not seem to affect  https://www.biostars.org/p/209118/ 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('T0570 vs T84:', dim(filtred_results)[1], sep=' ')), fill=TRUE)

sub = t2g[t2g$ensembl_gene_id %in% rownames(sorted_filtred_results), ]
sub = sub[, c(3, 4, 5, 7)]
subu <- sub[!duplicated(sub), ]
rownames(subu) <- subu$ensembl_gene_id

sorted_filtred_results$ensembl_id = subu[rownames(sorted_filtred_results), 'ensembl_gene_id']
sorted_filtred_results$gene_name = subu[rownames(sorted_filtred_results), 'external_gene_name']
sorted_filtred_results$biotype = subu[rownames(sorted_filtred_results), 'gene_biotype']
sorted_filtred_results$description2 = subu[rownames(sorted_filtred_results), 'description2']

# this is very dangerous, but lets go
sorted_filtred_results$median_clustera = rowMedians(subset[rownames(sorted_filtred_results), clustera_indexes]) 
sorted_filtred_results$median_clusterb = rowMedians(subset[rownames(sorted_filtred_results), clusterb_indexes])

write.table(sorted_filtred_results, file=paste(results_dir, '/effect_treatment_T2_0570.for.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, file=paste(results_dir, '/effect_treatment_T2_0570.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('genotype')) + ggtitle('effect treatment T2 0570')

#               
# volcano
#
plotting_x = sorted_filtred_results$log2FoldChange
y = sorted_filtred_results$padj
epsilon = min(y[y !=0])
plotting_y = -log10(y + epsilon) 
z = log10(rowMedians(txi$abundance[rownames(sorted_filtred_results), ]) + 1)
df = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20), plotting_z=clip(z, .min=0, .max=3))
reds = df[df$plotting_x > 0, ]
blues = df[df$plotting_x < 0, ]

plotting_x = anti_results$log2FoldChange
plotting_y = -log10(anti_results$padj)
blacks = data.frame(plotting_x=clip(plotting_x, .min=-6, .max=6), plotting_y=clip(plotting_y, .min=0, .max=20))

ggplot() + 
  geom_point(data=reds, aes(x=plotting_x, y=plotting_y, color=plotting_z), , size=3, shape=19, alpha=2/3, stroke=0) + 
  geom_point(data=blues, aes(plotting_x, plotting_y, color=plotting_z), size=3, shape=19, alpha=2/3, stroke=0) +
  geom_point(data=blacks, aes(plotting_x, plotting_y), color = "black", size=1, shape=19, alpha=0.2, stroke=0) +
  labs(x=expression('Expression [log'[2]~'FC]'), y=expression('Significance [log'[10]~'adjusted P]'), color=expression('Expression [log'[10]~'TPM]'), title='effect treatment T2 0570') + 
  theme_linedraw() +
  geom_segment(aes(x=-1, xend=-1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=1, xend=1, y=-log10(0.05), yend=20), linetype=2) +
  geom_segment(aes(x=-6, xend=-1, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  geom_segment(aes(x=1, xend=6, y=-log10(0.05), yend=-log10(0.05)), linetype=2) +
  xlim(-6, 6) +
  scale_color_viridis_c(option = "cividis") 
#ggsave(paste('effect_T0570_vs_T84', '.png', sep=''))
