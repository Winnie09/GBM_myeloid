rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
setwd("/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/")
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

# read in data
m_seu = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_subser.RDS')
expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M_log2cpm.rds')
saver <- readRDS('/home-4/zji4@jhu.edu/scratch/GBM/data/combine/saver/M.rds')
## select grade IV and untreated patients
meta =  m_seu@meta.data
saveRDS(meta, '/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/singleObject/M_meta.rds')
selectcell = rownames(meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ])
active.ident <- as.character(m_seu@active.ident)
names(active.ident) = rownames(meta)
saveRDS(active.ident, '/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/singleObject/M_active.ident.rds')
active.ident = active.ident[selectcell]

## read in diffusion mat which is reproduced on subsetted cell types
# phate = read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/Christina_share/subsetcelltypes/phate_sub.csv')

dm <- read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/data/GBM/Christina_share/m_subsetcelltypes/diffusion.csv', as.is = TRUE)

rownames(dm) = dm[,1]
dm = dm[,-1]
dm = dm[intersect(rownames(dm), names(active.ident)), ]

## pseudotime inference
library(TSCAN)
library(RColorBrewer)
# split up EMDSC cells into two cluster according to the distance to MMDSC and MAC1 s.t. the paths to MAC1 and MMDSC wouldn't have overlapping EMDSC cells.
# mac1center <- colMeans(dm[names(which(active.ident=='MAC1')),])
# mmdsccenter <- colMeans(dm[names(which(active.ident=='M-MDSC')),])
# mac1dist <- colSums((t(dm[names(which(active.ident=='E-MDSC')),])-mac1center)^2)
# mmdscdist <- colSums((t(dm[names(which(active.ident=='E-MDSC')),])-mmdsccenter)^2)
subident <- active.ident

# subident[names(which(mac1dist < mmdscdist))] <- 'E-MDSC_MAC1'
# subident[names(which(mac1dist > mmdscdist))] <- 'E-MDSC_MMDSC'

subident <- subident[rownames(dm)]

# pseudotime order
ordlist <- list()
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC','MAC1')]))
names(clu) = names(subident[subident %in% c('E-MDSC','MAC1')])
mc <- exprmclust(t(dm[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MAC1']] <- TSCANorder(mc,orderonly = T)
plotmclust(mc, cell_point_size = 0.1, x.lab = 'dm1', y.lab =  'dm2')
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC','MAC1','MAC2')]))
names(clu) = names(subident[subident %in% c('E-MDSC','MAC1','MAC2')])
mc <- exprmclust(t(dm[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MAC1_MAC2']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC','M-MDSC')]))
names(clu) = names(subident[subident %in% c('E-MDSC','M-MDSC')])
mc <- exprmclust(t(dm[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC','M-MDSC','PMN-MDSC')]))
names(clu) = names(subident[subident %in% c('E-MDSC','M-MDSC','PMN-MDSC')])
mc <- exprmclust(t(dm[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC_PMNMDSC']] <- TSCANorder(mc,orderonly = T)
#
clu = as.numeric(as.factor(subident[subident %in% c('E-MDSC','M-MDSC','MAC1')]))
names(clu) = names(subident[subident %in% c('E-MDSC','M-MDSC','MAC1')])
mc <- exprmclust(t(dm[names(clu),]),cluster=clu,reduce=F)
ordlist[['EMDSC_MMDSC_MAC1']] <- TSCANorder(mc,orderonly = T)
#
library(scattermore)
for (path in names(ordlist)){
  print(path)
  ord = ordlist[[path]]
  dir.create(paste0('./plot/', path), recursive = T)
  pdf(paste0('./plot/', path, '/pseudotime.pdf'), width = 5, height = 4)
  print(ggplot(data.frame(x = c(dm[ord,1], dm[!rownames(dm) %in% ord,1]), y = c(dm[ord,2], dm[!rownames(dm) %in% ord, 2]), time = c(seq(1, length(ord)), rep(NA, sum(!rownames(dm) %in% ord))))) + 
    geom_scattermore(aes(x = x, y = y, col = time), size = 0.05) + 
    scale_color_gradientn(colors = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(ord)), 'grey')) + 
    theme_classic() + xlab('DM1') + ylab('DM2'))
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
  
  expr.saver <- saver[rownames(expr.tmp), colnames(expr.tmp)]
  saveRDS(expr.saver, paste0('./result/', path, '/data/saver.rds'))
}


