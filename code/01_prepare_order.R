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
selectcell1 <- rownames(meta[meta$active.ident %in% c('E-MDSC', 'M-MDSC'), ])
selectcell2 <- rownames(meta[meta$active.ident %in% c('E-MDSC', 'MAC1'), ])


phate.tmp = phate[selectcell1,]
meta.tmp = meta[selectcell1,]
library(TSCAN)
library(RColorBrewer)
clu = as.character(meta.tmp$active.ident)
clu = as.numeric(as.factor(clu))
names(clu) = rownames(meta.tmp)
pd <- data.frame(x = phate.tmp[,1], y = phate.tmp[,2], clu = clu)
mc <- exprmclust(t(phate.tmp),cluster=clu,reduce=F)
plotmclust(mc,show_full_tree=T, cell_point_size = 0.01)
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
saveRDS(pseudotime, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MMDSC/pseudotime.rds')

expr.tmp <- as.matrix(expr[, selectcell1])
saveRDS(expr.tmp, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MMDSC/data/log2cpm.rds')

cellanno.tmp <- data.frame(cell = selectcell1, 
                           sample = meta[selectcell1, 'Patient'],
                           stringsAsFactors = FALSE)
saveRDS(cellanno.tmp, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MMDSC/cellanno.rds')



###########
phate.tmp = phate[selectcell2,]
meta.tmp = meta[selectcell2,]
expr.tmp <- as.matrix(expr[, selectcell2])
library(TSCAN)
library(RColorBrewer)
clu = as.character(meta.tmp$active.ident)
clu = as.numeric(as.factor(clu))
names(clu) = rownames(meta.tmp)
pd <- data.frame(x = phate.tmp[,1], y = phate.tmp[,2], clu = clu)
mc <- exprmclust(t(phate.tmp),cluster=clu,reduce=F)
plotmclust(mc,show_full_tree=T, cell_point_size = 0.01)
ord <- TSCANorder(mc, MSTorder=c(1,2), orderonly=T)
ggplot(data.frame(x = phate[ord,1], y = phate[ord,2], t = seq(1, length(ord)))) + 
  geom_point(aes(x = x, y = y, col = t), size = 0.1) + 
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, 'RdYlBu'))(length(ord))) + 
  theme_classic() + xlab('PHATE1') + ylab('PHATE2')   ## useful
pseudotime <- 1:length(ord)
names(pseudotime) <- ord
saveRDS(pseudotime, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MAC1/data/pseudotime.rds')
saveRDS(expr.tmp, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MAC1/data/log2cpm.rds')

cellanno.tmp <- data.frame(cell = selectcell2, 
                           sample = meta[selectcell2, 'Patient'],
                           stringsAsFactors = FALSE)
saveRDS(cellanno.tmp, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/EMDSC_MAC1/cellanno.rds')


