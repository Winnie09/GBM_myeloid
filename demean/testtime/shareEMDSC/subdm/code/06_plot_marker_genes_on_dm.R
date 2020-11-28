active.ident <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/singleObject/M_active.ident.rds')
dm <- read.csv('/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/Christina_share/m_subsetcelltypes/diffusion.csv', as.is = TRUE)
rownames(dm) = dm[,1]
dm = dm[,-1]

trajectory = 'EMDSC_MMDSC'
print(trajectory)
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/res/')
plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/plot/', trajectory, '/')

Res <- readRDS(paste0(rdir, 'ptest_res.rds'))
expr = Res$expr.demean
int = intersect(colnames(expr), rownames(dm))
dm = dm[int,]
expr = expr[, int]

library(ggplot2)
for (g in c('HK2', 'ENO2', 'MARCO', 'CD14', 'ENO2', 'CXCL2', 'S100A6', 'CXCL3', 'IL1B', 'VCAN')){
  pd = data.frame(DM1 = dm[,1], DM2 = dm[,2], centered.expression = expr[g,], stringsAsFactors = F)
  pdf(paste0(plotdir, 'dm_', g, '.pdf'), width = 4.2, height = 2.6)
  print(ggplot(data = pd, aes(x = DM1, y = DM2, color = centered.expression)) +
          geom_point(size = 0.2) +
          theme_classic() +
          scale_color_gradient2(low="blue",high="red", midpoint = 0) +
          ggtitle(g))
  dev.off()  
}

