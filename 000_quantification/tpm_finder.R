rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(stringr)

#
# 1. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/akthelia/results/000.quantification"
results_dir = kallisto_dir
# they are all a cases, --rf-stranded

#
# 1. generate gene to transcript mapping
#
df = read.csv('/Users/adrian/software/kallisto/human_index_standard/annotation.tsv', sep='\t')
t2g = df[, 2:3]
dim(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames = dirnames[grep('_processed', dirnames)]
paths = file.path(dirnames, 'kallisto_output_a/abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[8]))
labels = str_remove(labels, '_processed')

metadata = data.frame(labels)
metadata$path = paths
View(metadata)

#
# 3. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = metadata$labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, '/DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)
