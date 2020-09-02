rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
setwd("/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/")
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/treevar/subphate/plot/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/treevar/subphate/result/'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
# --------------------------------------------------------------
# input: seurat integrated object including:
# dimension reduction representation (e.g., umap, pca, phate)
# ct: dataframe/matrix, first column is cell name, second column is cell type, third column is sample.
# origin: the origin cell type
# --------------------------------------------------------------
# read in data
m_seu = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_subser.RDS')

## select grade IV and untreated patients
meta =  m_seu@meta.data
selectcell = rownames(meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ])
active.ident <- as.character(m_seu@active.ident)
names(active.ident) = rownames(meta)
phate = as.matrix(m_seu@reductions$phate@cell.embeddings)[selectcell, ]
active.ident = active.ident[selectcell]

## read in phate which is reproduced on subsetted cell types
phate = read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/Christina_share/subsetcelltypes/phate_sub.csv')
rownames(phate) = phate[,1]
phate = phate[,-1]
phate = phate[intersect(rownames(phate), names(active.ident)), ]

## tree variability
ct = data.frame(cell = rownames(phate), celltype = active.ident[rownames(phate)], sample = meta[rownames(phate), 'Patient'], stringsAsFactors = FALSE)
a = infer_tree_structure(pca = phate, ct = ct, origin.celltype = 'E-MDSC', plotdir = plotdir, number.cluster = length(unique(ct[,2]))+1, xlab = 'PHATE1', ylab = 'PHATE2')

result <- evaluate_uncertainty(a, 100)
saveRDS(result, paste0(rdir, 'result.rds'))

for (i in 1:length(result)){
  write.csv(result[[i]], paste0(rdir, names(result)[i], '.csv'), quote = FALSE, row.names = T, col.names = T)
}


