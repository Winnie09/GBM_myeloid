setwd('/scratch/users/whou10@jhu.edu/Wenpin/GBM_myeloid/demean')
library(Seurat)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
m_seu = readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_subser.RDS')
expr <- m_seu@assays$RNA@counts
rc <- colSums(expr)
rc <- rc/median(rc)
expr <- t(t(expr)/rc)
expr@x <- log2(expr@x + 1)
saveRDS(expr, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M_log2cpm.rds')

phate <- as.matrix(m_seu@reductions$phate@cell.embeddings)
umap <- as.matrix(m_seu@reductions$umap@cell.embeddings)
pca <- as.matrix(m_seu@reductions$pca@cell.embeddings)
meta <- m_seu@meta.data
active.ident <- as.character(m_seu@active.ident)
library(ggplot2)
library(scattermore)
ggplot(data = data.frame(x = phate[,1], y = phate[,2], celltype = active.ident), 
       aes(x = x, y = y, color = active.ident)) + 
  geom_scattermore()

meta = cbind(meta, active.ident = active.ident)
meta <- meta[meta$Tumor.Grade == 'IV' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('E-MDSC', 'M-MDSC', 'MAC1'), ])

phate.tmp = phate[selectcell,]
meta.tmp = meta[selectcell,]
library(TSCAN)
library(RColorBrewer)
clu = as.character(meta.tmp$active.ident)
clu = as.numeric(as.factor(clu))
names(clu) = rownames(meta.tmp)
pd <- data.frame(x = phate.tmp[,1], y = phate.tmp[,2], clu = clu)
mc <- exprmclust(t(phate.tmp),cluster=clu,reduce=F)
# plotmclust(mc,show_full_tree=T, cell_point_size = 0.01)

## EMDSC -> MMDSC
ord <- TSCANorder(mc, MSTorder=c(1,2), orderonly=T)
ggplot(data.frame(ct=meta[ord, ]$active.ident,
                  pt=1:length(ord)),
       aes(x=pt,y=ct)) + 
  geom_point(size = 0.1)
ggplot(data.frame(x = phate[ord,1], y = phate[ord,2], t = seq(1, length(ord)))) + 
  geom_point(aes(x = x, y = y, col = t), size = 0.1) + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord))) + 
  theme_classic() + xlab('PHATE1') + ylab('PHATE2')   ## useful
pseudotime <- 1:length(ord)
names(pseudotime) <- ord
dir.create('./result/EMDSC_MMDSC/data', recursive = T)
saveRDS(pseudotime, './result/EMDSC_MMDSC/data/pseudotime.rds')

expr.tmp <- as.matrix(expr[, ord])
saveRDS(expr.tmp, './result/EMDSC_MMDSC/data/log2cpm.rds')

cellanno.tmp <- data.frame(cell = ord, 
                           sample = meta[ord, 'Patient'],
                           stringsAsFactors = FALSE)
saveRDS(cellanno.tmp, './result/EMDSC_MMDSC/data/cellanno.rds')


## EMDSC -> MAC1
ord2 <- TSCANorder(mc, MSTorder=c(1,3), orderonly=T)
ggplot(data.frame(x = phate[ord2,1], y = phate[ord2,2], t = seq(1, length(ord2)))) + 
  geom_point(aes(x = x, y = y, col = t), size = 0.1) + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord2))) + 
  theme_classic() + xlab('PHATE1') + ylab('PHATE2')   ## useful
pseudotime <- 1:length(ord2)
names(pseudotime) <- ord2
dir.create('./result/EMDSC_MAC1/data', recursive = T)
saveRDS(pseudotime, './result/EMDSC_MAC1/data/pseudotime.rds')

expr.tmp <- as.matrix(expr[, ord2])
saveRDS(expr.tmp, './result/EMDSC_MAC1/data/log2cpm.rds')
cellanno.tmp <- data.frame(cell = ord2, 
                           sample = meta[ord2, 'Patient'],
                           stringsAsFactors = FALSE)
saveRDS(cellanno.tmp, './result/EMDSC_MAC1/data/cellanno.rds')

