rm(list = ls())
library(ggplot2)
library(reshape2)
library(RColorBrewer)
suppressMessages(library(igraph))
library(parallel)
library(splines)
library(viridis)
library(here)

source(
  '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R'
)
trajectory = 'EMDSC_MMDSC'
trajectory = as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
# for (trajectory in c('EMDSC_MAC1','EMDSC_MAC1_MAC2','EMDSC_MMDSC','EMDSC_MMDSC_PMNMDSC','EMDSC_MMDSC_MAC1')) {
print(trajectory)
rdir <-
  here('notdemean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'result',
       trajectory,
       'res')
plotdir <-
  here('notdemean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'plot',
       trajectory)
dir.create(plotdir, recursive = T)

Res <- readRDS(paste0(rdir, '/ptest_res.rds'))
res = data.frame(
  fdr = Res$fdr,
  foldchange = Res$foldchange,
  pvalue = Res$pvalue
)
res = res[res[, 1] < 0.05, ]
res = res[order(res[, 1], res[, 2]), ]
str(res)
write.csv(res, paste0(rdir, '/differential_genes.csv'))

pt = Res$pseudotime
cellanno = Res$cellanno
active.ident <-
  readRDS(
    '/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/singleObject/Myeloid/active.ident.rds'
  )

## -----------------
## cluster genes
## -----------------
# [option 1/2] use population pattern for genes: phi * x * beta
pred <- t(sapply(rownames(res), function(i) {
  get_population_fit(testobj = Res,
                     variable = NA,
                     gene = i)
}))
## [option 2/2] use demean fitted values for genes
# pd = Res$expr.demean[rownames(res),]
# pred  <- t(sapply(1:nrow(pd), function(i) {
#   loess(pd[i, ] ~ seq(1, ncol(pd)))$fitted
# }))
#
# dimnames(pred) <- dimnames(Res$expr.ori[rownames(res),])
saveRDS(pred, paste0(rdir, '/pred.rds'))


## standardize the values
cm <- rowMeans(pred)
csd <- apply(pred, 1, sd)
pred.scale <- (pred - cm) / csd
set.seed(12345)
if (trajectory == 'EMDSC_MMDSC_PMNMDSC') {
  clu <- kmeans(pred.scale, 7)$cluster
} else {
  clu <- kmeans(pred.scale, 5)$cluster
}
clu <- sort(clu)
res$clu = clu[rownames(res)]
res$cor <- sapply(rownames(res), function(i) cor(pred.scale[i,], seq(1,ncol(pred.scale))))
pred.scale <- pred.scale[rownames(res)[order(res$clu, res$cor)], ]
write.csv(res, paste0(rdir, '/differential_genes_clu.csv'))

## -------
## plot
## -------
# annotate rows and columns
colann <-
  data.frame(sample = cellanno[match(cellanno[, 1], colnames(pred.scale)), 2],
             celltype = active.ident[colnames(pred.scale)],
             stringsAsFactors = F)
rownames(colann) = colnames(pred.scale)
rowann = data.frame(cluster = as.character(clu), stringsAsFactors = F)
rownames(rowann) = names(clu)

# define colors
library(pheatmap)
col.clu = brewer.pal(length(unique(clu)), 'Set1')
names(col.clu) = unique(clu)
col.ct = brewer.pal(length(unique(colann$celltype)), name = 'Dark2')
names(col.ct) = unique(colann$celltype)
col.ct = col.ct[!is.na(names(col.ct))]

col.sample = colorRampPalette(rev(brewer.pal(n = 7, name = "Set1")))(length(unique(colann$sample)))
names(col.sample) = unique(colann$sample)

if (trajectory == 'EMDSC_MMDSC_PMNMDSC') {
  png(
    paste0(plotdir, '/hm_kmeans_scale.png'),
    width = 1575,
    height = 2400,
    res = 250
  )
} else {
  png(
    paste0(plotdir, '/hm_kmeans_scale.png'),
    width = 1050,
    height = 1600,
    res = 200
  )
}
cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MAC1
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MAC1_MAC2
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 60)) ## EMDSC_MMDSC
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MMDSC_PMNMDSC
# cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MMDSC-MAC1
pheatmap(
  pred.scale,
  cluster_rows = F,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = cpl,
  annotation_row = rowann,
  annotation_col = colann,
  annotation_colors = list(
    cluster = col.clu,
    celltype = col.ct,
    sample = col.sample
  ),
  cellwidth = 241.24 / ncol(pred.scale),
  cellheight = 467.35 / nrow(pred.scale)
)
dev.off()


agg <-
  t(sapply(unique(clu), function(i)
    colMeans(pred.scale[clu == i, ])))
rownames(agg) <- paste0('cluster', unique(clu))
library(ggplot2)
pd <- reshape2::melt(agg)
pd$pseudotime = pt[pd[, 2]]
colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
pd$cluster = paste0(as.character(pd$cluster), '(', table(clu)[gsub('cluster','',as.character(pd$cluster))],')')
cluord <- order(as.numeric(sapply(unique(pd$cluster), function(i) sub('\\(.*', '',sub('cluster', '', i)))))
pd$cluster <-
  factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])


pdf(paste0(plotdir, '/gene_cluster_pattern.pdf'),
    width = 4.2,
    height = 3)
print(
  ggplot(
    data = pd,
    aes(
      x = pseudotime,
      y = expression,
      group = cluster,
      color = cluster
    )
  ) +
    geom_smooth() +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    xlab('Pseudotime') +
    ylab('Averaged Fitted Expression')
)
dev.off()


gl <- c(
  'CD14',
  'CD44',
  'CD52',
  'CXCL5',
  'CD36',
  'IL1B',
  'CCL7',
  'CCL4',
  'FTH1',
  'HLA-DQB1',
  'CXCL2',
  'CXCL3',
  'ENO2',
  'HK2',
  'MARCO'
)
fit = pred.scale[rownames(res)[abs(res$cor)>0.3],]
colnames(fit) <- pt[colnames(pred.scale)]
png(paste0(plotdir, '/smooth_hm.png'),
    width = 800,
    height = 1000)
mySTIP(fit,
     gl)
dev.off()


pdf(paste0(plotdir, '/smooth_hm.pdf'),
    width = 8,
    height = 10)
mySTIP(fit,
     gl)
dev.off()

## plot marker gene on dm
dm <- read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/Christina_share/m_subsetcelltypes/diffusion.csv', as.is = TRUE)
rownames(dm) = dm[,1]
dm = dm[,-1]

int = intersect(colnames(pred), rownames(dm))
dm = dm[int,]

for (g in gl){
  pd = data.frame(DM1 = dm[,1], DM2 = dm[,2], Expression = pred[g,], stringsAsFactors = F)
  pdf(paste0(plotdir, '/dm_', g, '.pdf'), width = 3.5, height = 2.6)
  print(ggplot(data = pd, aes(x = DM1, y = DM2, color = Expression)) +
          geom_point(size = 0.1, alpha = 0.5) +
          theme_classic() +
          scale_color_gradientn(colors = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(nrow(dm)), 'grey')) + 
          ggtitle(g))
  dev.off()  
}

  
