library(ggplot2)
library(gridExtra)
library(here)
here()
for (traj in list.files(here('demean', 'testtime', 'seperateEMDSC', 'subdm', 'result'))) {
  print(traj)
  a = read.csv(
    here('demean', 'testtime', 'seperateEMDSC', 'subphate', 'result', traj, 'res', 'differential_genes_clu.csv'),
    as.is = T
  )
  
  b = read.csv(
    here('demean', 'testtime', 'seperateEMDSC', 'subdm', 'result', traj, 'res', 'differential_genes_clu.csv'),
    as.is = T
  )
  
  c = read.csv(
    here('demean', 'testtime', 'shareEMDSC','subdm', 'result', traj, 'res', 'differential_genes_clu.csv'),
    as.is = T
  )
  
  str(a)
  str(b)
  str(c)
  
  length(intersect(a[, 1], b[, 1]))
  
  length(intersect(a[, 1], c[, 1]))
  
  length(intersect(c[, 1], b[, 1]))
  
  t1 <-
    readRDS(
      here('demean','testtime','seperateEMDSC','subphate','result',traj,'data','pseudotime.rds')
    )
  
  t2 <-
    readRDS(
      here('demean','testtime','seperateEMDSC','subdm','result',traj,'data','pseudotime.rds')
    )
  
  t3 <-
    readRDS(
      here('demean','testtime','shareEMDSC','subdm','result',traj,'data','pseudotime.rds')
    )

  str(t1)
  str(t2)
  str(t3)  
  int = union(union(names(t1), names(t2)), names(t3))
  uni12 = intersect(names(t1), names(t2))
  uni13 = intersect(names(t1), names(t3))
  uni23 = intersect(names(t2), names(t3))
  str(uni12)
  str(uni23)
  str(uni13)
  
  pd <-
    data.frame(
      t1 = t1[int],
      t2 = t2[int],
      t3 = t3[int],
      stringsAsFactors = F
    )
  
  p1 <- ggplot(data = pd) +
    geom_point(aes(x = t1, y = t2), size = 0.1, alpha = 0.2) +
    ggtitle(paste0('PCC=',  round(cor(t1[uni12], t2[uni12]), 3))) +
    xlab('seperateEMDSC/subphate') + ylab('seperateEMDSC/subdm')
  
  p2 <- ggplot(data = pd) +
    geom_point(aes(x = t1, y = t3), size = 0.1, alpha = 0.2) +
    ggtitle(paste0('PCC=',  round(cor(t1[uni13], t1[uni13]), 3))) +
    xlab('seperateEMDSC/subphate') + ylab('shareEMDSC/subdm')
  
  p3 <- ggplot(data = pd) +
    geom_point(aes(x = t2, y = t3), size = 0.1, alpha = 0.2) +
    ggtitle(paste0('PCC=',  round(cor(t2[uni23], t3[uni23]), 3))) +
    xlab('seperateEMDSC/subdm') + ylab('shareEMDSC/subdm')
  
  pdf(
    here(
      'demean',
      'testtime',
      'diagnosis',
      'plot',
      paste0('compare_pseudotime_', traj, '.pdf')
    ),
    width = 8,
    height = 2.5
  )
  print(grid.arrange(p1, p2, p3, nrow = 1))
  dev.off()
}

