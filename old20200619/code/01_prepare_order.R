suppressMessages(library(Seurat))
obj <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/seuratObject/M_ser.RDS')
cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/count.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')

mdscProp <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv')
active.ident = obj@active.ident
pca <- obj@reductions$pca@cell.embeddings
umap <- obj@reductions$umap@cell.embeddings
phate <- obj@reductions$phate@cell.embeddings
saveRDS(active.ident, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/active.ident.rds')
saveRDS(pca, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/pca.rds')
saveRDS(umap, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/umap.rds')
saveRDS(phate, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/phate.rds')

DimPlot(obj,label=T)


meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Grade == 'High' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('MDSC', 'MAC1', 'MAC2', 'MAC3', 'NEU1', 'MAC4'), ])

library(ggplot2)
ggplot() + geom_point(data = data.frame(x = pca[selectcell, 1], y = pca[selectcell, 2], ct = as.factor(active.ident[selectcell])), aes(x = x, y = y, col=ct), size=0.01) 

pca <- pca[, 1:10]
library(umap)
set.seed(13)
u <- umap(d = pca[selectcell, 1:10])$layout
saveRDS(u, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/umap.rds')
ggplot() + geom_point(data = data.frame(x = u[,1], y = u[,2], ct = as.factor(active.ident[rownames(u)])), aes(x = x, y = y, col = ct), size = 0.01) + xlab('UMAP1') + ylab('UMAP2') + guides(color=guide_legend(override.aes = list(size=5)))

### try to recover the old UMAP ==
for (s in 5:100){
  print(s)
  set.seed(s)
  u <- umap(d = pca[selectcell, 1:10])$layout
  pdf(paste0(s,'.pdf'), width = 5, height = 5)
  print(ggplot() + geom_point(data = data.frame(x = u[,1], y = u[,2], ct = as.factor(active.ident[rownames(u)])), aes(x = x, y = y, col = ct), size = 0.01) + xlab('UMAP1') + ylab('UMAP2') + guides(color=guide_legend(override.aes = list(size=5))))
  dev.off()
}
###  
library(TSCAN)
set.seed(12345)
clu <- kmeans(pca[selectcell,], 5)$cluster
ggplot() + geom_point(data = data.frame(x = u[,1], y = u[,2], clu = as.factor(clu[rownames(u)])), aes(x = x, y = y, col = clu), size = 0.01) + xlab('UMAP1') + ylab('UMAP2') + guides(color=guide_legend(override.aes = list(size=5))) + theme_classic()

tmpclu <- clu[clu %in% c(5, 3, 1)]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <-  n
order <- rev(TSCANorder(exprmclust(t(pca[n,]),cluster=tmpclu,reduce=F),orderonly = T))
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_NEU1.rds')

tmpclu <- clu[clu %in% c(5, 4)]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <-  n
order <- rev(TSCANorder(exprmclust(t(pca[n,]),cluster=tmpclu,reduce=F),orderonly = T))
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC1.rds')

tmpclu <- clu[clu %in% c(2, 4)]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <-  n
order <- TSCANorder(exprmclust(t(pca[n,]),cluster=tmpclu,reduce=F),orderonly = T)
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MAC2_MAC1.rds')

clu <- as.character(active.ident[selectcell])
tmpclu <- clu[clu %in% c('MDSC', 'MAC3','MAC1')]
names(tmpclu) <-  selectcell[clu %in% c('MDSC', 'MAC3','MAC1')]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <- n
set.seed(123)
order <- rev(TSCANorder(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F),orderonly = T, MSTorder=c(1,2,3)))
# plotmclust(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F))
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_MAC1.rds')


clu <- as.character(active.ident[selectcell])
tmpclu <- clu[clu %in% c('MDSC', 'MAC3')]
names(tmpclu) <-  selectcell[clu %in% c('MDSC', 'MAC3')]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <- n
set.seed(123)
order <- rev(TSCANorder(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F),orderonly = T, MSTorder=c(1,2)))
# plotmclust(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F))
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3.rds')

###

tmpclu <- clu[clu %in% c( 'MAC3', 'MAC1')]
names(tmpclu) <-  selectcell[clu %in% c('MAC3', 'MAC1')]
n <- names(tmpclu)
tmpclu <- as.numeric(as.factor(tmpclu))
names(tmpclu) <- n
set.seed(123)
order <- rev(TSCANorder(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F),orderonly = T, MSTorder=c(1,2)))
# plotmclust(exprmclust(t(pca[names(tmpclu),]),cluster=tmpclu,reduce=F))
ggplot() + geom_point(data=data.frame(x = u[order,1], y = u[order, 2], pseudotime = seq(1, length(n))), aes(x = x, y = y, col=pseudotime), size = 0.1, alpha = 0.5) + scale_color_gradient(low = 'grey', high = 'red') + xlab('UMAP1') + ylab('UMAP2') + theme_classic()
saveRDS(order, '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MAC3_MAC1.rds')

