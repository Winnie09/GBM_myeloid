rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
setwd("/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subphate/")
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

# read in data
m_seu = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_subser.RDS')
expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M_log2cpm.rds')

## select grade IV and untreated patients
meta =  m_seu@meta.data
saveRDS(meta, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M_meta.rds')
selectcell = rownames(meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ])
active.ident <- as.character(m_seu@active.ident)
names(active.ident) = rownames(meta)
saveRDS(active.ident, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M_active.ident.rds')
active.ident = active.ident[selectcell]

## read in phate which is reproduced on subsetted cell types
phate = read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/Christina_share/subsetcelltypes/phate_sub.csv')
rownames(phate) = phate[,1]
phate = phate[,-1]
phate = phate[intersect(rownames(phate), names(active.ident)), ]

## pseudotime inference
library(TSCAN)
library(RColorBrewer)
# split up EMDSC cells into two cluster according to the distance to MMDSC and MAC1 s.t. the paths to MAC1 and MMDSC wouldn't have overlapping EMDSC cells.
mac1center <- colMeans(phate[names(which(active.ident=='MAC1')),])
mmdsccenter <- colMeans(phate[names(which(active.ident=='M-MDSC')),])
mac1dist <- colSums((t(phate[names(which(active.ident=='E-MDSC')),])-mac1center)^2)
mmdscdist <- colSums((t(phate[names(which(active.ident=='E-MDSC')),])-mmdsccenter)^2)
subident <- active.ident

subident[names(which(mac1dist < mmdscdist))] <- 'E-MDSC_MAC1'
subident[names(which(mac1dist > mmdscdist))] <- 'E-MDSC_MMDSC'

subident <- subident[rownames(phate)]

# pseudotime order
ordlist <- list()
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC_MAC1','MAC1')]))
names(clu) = names(subident[subident %in% c('E-MDSC_MAC1','MAC1')])
mc <- exprmclust(t(phate[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MAC1']] <- TSCANorder(mc,orderonly = T)
plotmclust(mc, cell_point_size = 0.1, x.lab = 'PHATE1', y.lab =  'PHATE2')
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC_MAC1','MAC1','MAC2')]))
names(clu) = names(subident[subident %in% c('E-MDSC_MAC1','MAC1','MAC2')])
mc <- exprmclust(t(phate[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MAC1_MAC2']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC')]))
names(clu) = names(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC')])
mc <- exprmclust(t(phate[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC','PMN-MDSC')]))
names(clu) = names(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC','PMN-MDSC')])
mc <- exprmclust(t(phate[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC_PMNMDSC']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC','MAC1')]))
names(clu) = names(subident[subident %in% c('E-MDSC_MMDSC','M-MDSC','MAC1')])
mc <- exprmclust(t(phate[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC-MAC1']] <- TSCANorder(mc,orderonly = T)
#

for (path in names(ordlist)){
  print(path)
  ord = ordlist[[path]]
  dir.create(paste0('./plot/', path), recursive = T)
  pdf(paste0('./plot/', path, '/pseudotime.pdf'), width = 5, height = 4)
  print(ggplot(data.frame(x = c(phate[ord,1], phate[!rownames(phate) %in% ord,1]), y = c(phate[ord,2], phate[!rownames(phate) %in% ord, 2]), time = c(seq(1, length(ord)), rep(NA, sum(!rownames(phate) %in% ord))))) + 
    geom_point(aes(x = x, y = y, col = time), size = 0.05) + 
    scale_color_gradientn(colors = c(colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord)), 'grey')) + 
    theme_classic() + xlab('PHATE1') + ylab('PHATE2'))
  dev.off()
  pseudotime <- 1:length(ord)
  names(pseudotime) <- ord
  dir.create(paste0('./result/', path, '/data'), recursive = T)
  saveRDS(pseudotime, paste0('./result/', path, '/data/pseudotime.rds'))
  
  expr.tmp <- as.matrix(expr[, ord])
  saveRDS(expr.tmp, paste0('./result/', path, '/data/log2cpm.rds'))

  cellanno.tmp <- data.frame(cell = ord, 
                           sample = meta[ord, 'Patient'],
                           stringsAsFactors = FALSE)
  saveRDS(cellanno.tmp, paste0('./result/', path, '/data/cellanno.rds'))
}


