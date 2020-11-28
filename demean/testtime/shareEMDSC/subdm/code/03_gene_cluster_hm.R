rm(list=ls())
library(ggplot2)
library(Seurat)
library(reshape2)
library(TSCAN)
library(scattermore)
library(RColorBrewer)
suppressMessages(library(igraph))
library(parallel)
library(splines)
library(viridis)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

for (trajectory in c('EMDSC_MAC1', 'EMDSC_MAC1_MAC2', 'EMDSC_MMDSC', 'EMDSC_MMDSC_PMNMDSC', 'EMDSC_MMDSC_MAC1')){
  print(trajectory)
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/res/')
  plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/plot/', trajectory, '/')
  
  Res <- readRDS(paste0(rdir, 'ptest_res.rds'))
  pt = Res$pseudotime
  
  
  res <- read.csv(paste0(rdir, 'differential_genes_clu.csv'), row.names = 1)
  clu = res$clu
  names(clu) = rownames(res)
  
  pred <- Res$predict.values[rownames(res), names(pt)]
  agg <- t(sapply(unique(clu), function(i) colMeans(pred[clu == i, ,drop=F])))
  rownames(agg) <- paste0('cluster', unique(clu), '(',table(clu), ')')
  
  library(ggplot2)
  pd <- reshape2::melt(agg)
  pd$pseudotime = pt[pd[,2]]
  colnames(pd) <- c('cluster', 'cell', 'expression', 'pseudotime')
  
  cluord <- order(as.numeric(gsub('\\(.*', '',sub('cluster', '', unique(pd$cluster)))))
  pd$cluster <- factor(as.character(pd$cluster), levels =  unique(pd$cluster)[cluord])
    
  pdf(paste0(plotdir, 'gene_cluster_pattern.pdf'), width = 4.5, height = 3.2)  
  print(ggplot(data= pd, aes(x = pseudotime, y = expression, group = cluster, color = cluster)) + 
    geom_smooth() +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    xlab('Pseudotime') +
    ylab('Averaged Fitted Centered Expression'))
  dev.off()
}
  

