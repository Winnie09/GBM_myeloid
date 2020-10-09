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

#trajectory <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
for (trajectory in c('EMDSC_MAC1','EMDSC_MAC1_MAC2','EMDSC_MMDSC','EMDSC_MMDSC_PMNMDSC','EMDSC_MMDSC-MAC1')) {
  print(trajectory)
  rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/res/')
  datadir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/data/')
  plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/plot/', trajectory, '/')
  
  Res <- readRDS(paste0(rdir, 'ptest_res.rds'))
  res = data.frame(fdr = Res$fdr, foldchange = Res$foldchange, pvalue = Res$pvalue)
  res = res[res[,1]<0.05, ]
  res = res[order(res[,1], res[,2]), ]
  str(res)
  write.csv(res, paste0(rdir, 'differential_genes.csv'))
  
  pt = Res$pseudotime
  pred <- Res$predict.values[rownames(res), names(pt)]
  str(pred)
  cellanno = Res$cellanno
  
  active.ident <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/singleObject/M_active.ident.rds') ###
  
  
  # cluster genes
  set.seed(12345)
  clu <- kmeans(pred,5)$cluster
  clu <- sort(clu)
  table(clu)
  res$clu = clu[rownames(res)]
  write.csv(res, paste0(rdir, 'differential_genes_clu.csv'))
  
  # annotate rows and columns
  colann <- data.frame(sample = cellanno[match(cellanno[,1], colnames(pred)),2],
                       celltype = active.ident[colnames(pred)],stringsAsFactors = F)
  rownames(colann) = colnames(pred)
  rowann = data.frame(cluster = as.character(clu),stringsAsFactors = F)
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
  
  if (trajectory == 'EMDSC_MMDSC_PMNMDSC'){
    png(paste0(plotdir, 'hm_kmeans.png'),width=1575,height=2400, res=250)
  } else {
    png(paste0(plotdir, 'hm_kmeans.png'),width=1050,height=1600, res=200)
  }
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  # cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MAC1
  # cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MAC1_MAC2
  # cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 60)) ## EMDSC_MMDSC
  # cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MMDSC_PMNMDSC
  cpl = c(rep(cpl[1], 100), cpl, rep(cpl[length(cpl)], 100)) ## EMDSC_MMDSC-MAC1
  
  pheatmap(pred[names(clu), ], 
           cluster_rows = F, cluster_cols = FALSE,
           show_rownames = FALSE, show_colnames = FALSE,
           color = cpl,
           annotation_row = rowann, annotation_col = colann,
           annotation_colors = list(cluster = col.clu, celltype = col.ct, sample = col.sample),
           cellwidth = 241.24/ncol(pred), cellheight = 467.35/nrow(pred))
  dev.off()
}


