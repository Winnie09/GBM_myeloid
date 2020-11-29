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
trajectory = 'EMDSC_MAC1'
# for (trajectory in c('EMDSC_MAC1','EMDSC_MAC1_MAC2','EMDSC_MMDSC','EMDSC_MMDSC_PMNMDSC','EMDSC_MMDSC_MAC1')) {
print(trajectory)
rdir <-
  here('demean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'result',
       trajectory,
       'res')
datadir <-
  here('demean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'result',
       trajectory,
       'data')
plotdir <-
  here('demean', 'testtime', 'shareEMDSC', 'subdm', 'plot', trajectory)

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

# cluster genes
## use population pattern for genes: phi * x * beta
# pred <- t(sapply(rownames(res), function(i){
#   get_population_fit(testobj = Res, variable = NA, gene = i)
# }))
pd = Res$expr.demean[rownames(res),]
fit.demean <- t(sapply(1:nrow(pd), function(i) {
  loess(pd[i, ] ~ seq(1, ncol(pd)))$fitted
}))

pred <- fit.demean
dimnames(pred) <- dimnames(Res$expr.demean[rownames(res),])
## standardize the values
cm <- rowMeans(pred)
csd <- apply(pred, 1, sd)
pred <- (pred - cm) / csd
set.seed(12345)
clu <- kmeans(pred, 5)$cluster
clu <- sort(clu)
res$clu = clu[rownames(res)]
res$cor <- sapply(rownames(res), function(i) cor(pred[i,], seq(1,ncol(pred))))
pred <- pred[rownames(res)[order(res$clu, res$cor)], ]

# write.csv(res, paste0(rdir, 'differential_genes_clu.csv'))

# annotate rows and columns
colann <-
  data.frame(sample = cellanno[match(cellanno[, 1], colnames(pred)), 2],
             celltype = active.ident[colnames(pred)],
             stringsAsFactors = F)
rownames(colann) = colnames(pred)
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
    paste0(plotdir, 'hm_kmeans_scale.png'),
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
  pred,
  cluster_rows = F,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE
)

pheatmap(
  pred,
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
  cellwidth = 241.24 / ncol(pred),
  cellheight = 467.35 / nrow(pred)
)
dev.off()
# }



agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ])))
  rownames(agg) <- paste0('cluster', unique(clu))
  
  library(ggplot2)
  pd <- reshape2::melt(agg)
  pd$pseudotime = pt[pd[,2]]
  colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
  pd$cluster <- factor(as.character(pd$cluster), levels =  paste0('cluster', 1:5))
    
  print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
    geom_smooth() +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    xlab('Pseudotime') +
    ylab('Averaged Fitted Centered Expression'))
  dev.off()


