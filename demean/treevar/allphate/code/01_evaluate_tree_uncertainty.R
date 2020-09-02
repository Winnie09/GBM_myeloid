rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
setwd("/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/treevar/allphate/")
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/treevar/allphate/plot/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/treevar/allphate/result/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
# --------------------------------------------------------------
# input: seurat integrated object including:
# umap, pca
# ct: dataframe/matrix, first column is cell name, second column is cell type, third column is sample.
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
m_seu = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_subser.RDS')
meta =  m_seu@meta.data
selectcell = rownames(meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ])

active.ident <- as.character(m_seu@active.ident)
names(active.ident) = rownames(meta)

phate = as.matrix(m_seu@reductions$phate@cell.embeddings)[selectcell, ]
active.ident = active.ident[selectcell]

set.seed(12345)
id <- sample(1:nrow(phate), 1e4)
table(active.ident[id])
active.ident = active.ident[id]
phate = phate[id,]
ct = data.frame(cell = names(active.ident), celltype = active.ident, sample = meta[names(active.ident), 'Patient'], stringsAsFactors = FALSE)

a = infer_tree_structure(pca = phate, ct = ct, origin.celltype = 'E-MDSC', plotdir = plotdir, number.cluster = length(unique(active.ident)), xlab = 'PHATE1', ylab = 'PHATE2')
result <- evaluate_uncertainty(a, 100)
saveRDS(result, paste0(rdir, 'result.rds'))
for (i in 1:length(result)){
  write.csv(result[[i]], paste0(rdir, names(result)[i], '.csv'), quote = FALSE, row.names = T, col.names = T)
}


