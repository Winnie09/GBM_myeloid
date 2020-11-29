library(ggplot2)
library(gridExtra)
library(here)
here()
# for (traj in list.files(here('demean', 'testtime', 'seperateEMDSC', 'subdm', 'result'))) {
traj = 'EMDSC_MMDSC'
print(traj)
a = read.csv(
  here(
    'demean',
    'testtime',
    'seperateEMDSC',
    'subphate',
    'result',
    traj,
    'res',
    'differential_genes_clu.csv'
  ),
  as.is = T
)

b = read.csv(
  here(
    'demean',
    'testtime',
    'seperateEMDSC',
    'subdm',
    'result',
    traj,
    'res',
    'differential_genes_clu.csv'
  ),
  as.is = T
)

c = read.csv(
  here(
    'demean',
    'testtime',
    'shareEMDSC',
    'subdm',
    'result',
    traj,
    'res',
    'differential_genes_clu.csv'
  ),
  as.is = T
)


for (clu in 1:5) {
  print(paste0('cluster', clu))
  print(length(intersect(a[a[, 5] == 5, 1], b[b[, 5] == clu, 1])))
  print('\n')
  print(length(intersect(a[a[, 5] == 5, 1], c[c[, 5] == clu, 1])))
}

selgene <- a[a[, 5] == 5, 1]
g = selgene[1]
a_res <- readRDS(
  here(
    'demean',
    'testtime',
    'seperateEMDSC',
    'subphate',
    'result',
    traj,
    'res',
    'ptest_res.rds'
  )
)
names(a_res)

for (g in selgene) {
  pd <-
    data.frame(
      Expression = a_res$expression[g, names(a_res$pseudotime)],
      Pseudotime = a_res$pseudotime,
      stringsAsFactors = F
    )
  pdf(
    here(
      'demean',
      'testtime',
      'seperateEMDSC',
      'subphate',
      'plot',
      traj,
      paste0(g, '.pdf')
    ),
    width = 3,
    height = 3
  )
  print(ggplot(pd, aes(x = Pseudotime, y = Expression)) +
          geom_point(size = 0.1, alpha = 0.2) +
          geom_smooth() +
          theme_classic())
  dev.off()
}

b_res <- readRDS(
  here(
    'demean',
    'testtime',
    'seperateEMDSC',
    'subdm',
    'result',
    traj,
    'res',
    'ptest_res.rds'
  )
)
names(b_res)
for (g in selgene) {
  pd <-
    data.frame(
      Expression = b_res$expr.demean[g, names(b_res$pseudotime)],
      Pseudotime = b_res$pseudotime,
      stringsAsFactors = F
    )
  pdf(
    here(
      'demean',
      'testtime',
      'seperateEMDSC',
      'subdm',
      'plot',
      traj,
      paste0(g, '.pdf')
    ),
    width = 3,
    height = 3
  )
  print(ggplot(pd, aes(x = Pseudotime, y = Expression)) +
          geom_point(size = 0.1, alpha = 0.2) +
          geom_smooth() +
          theme_classic())
  dev.off()
  
}


c_res <- readRDS(
  here(
    'demean',
    'testtime',
    'shareEMDSC',
    'subdm',
    'result',
    traj,
    'res',
    'ptest_res.rds'
  )
)
names(c_res)

for (g in selgene) {
  pd <-
    data.frame(
      Expression = c_res$expr.demean[g, names(c_res$pseudotime)],
      Pseudotime = c_res$pseudotime,
      stringsAsFactors = F
    )
  pdf(
    here(
      'demean',
      'testtime',
      'shareEMDSC',
      'subdm',
      'plot',
      traj,
      paste0(g, '_demean.pdf')
    ),
    width = 3,
    height = 3
  )
  print(ggplot(pd, aes(x = Pseudotime, y = Expression)) +
          geom_point(size = 0.1, alpha = 0.2) +
          geom_smooth() +
          theme_classic())
  dev.off()
  
  pd <-
    data.frame(
      Expression = c_res$expr.ori[g, names(c_res$pseudotime)],
      Pseudotime = c_res$pseudotime,
      stringsAsFactors = F
    )
  pdf(
    here(
      'demean',
      'testtime',
      'shareEMDSC',
      'subdm',
      'plot',
      traj,
      paste0(g, '_ori.pdf')
    ),
    width = 3,
    height = 3
  )
  print(ggplot(pd, aes(x = Pseudotime, y = Expression)) +
          geom_point(size = 0.1, alpha = 0.2) +
          geom_smooth() +
          theme_classic())
  dev.off()

  
  pd <-
    data.frame(
      Expression = c_res$predict.values[g, names(c_res$pseudotime)],
      Pseudotime = c_res$pseudotime,
      stringsAsFactors = F
    )
  pdf(
    here(
      'demean',
      'testtime',
      'shareEMDSC',
      'subdm',
      'plot',
      traj,
      paste0(g, '_pred.pdf')
    ),
    width = 3,
    height = 3
  )
  print(ggplot(pd, aes(x = Pseudotime, y = Expression)) +
          geom_point(size = 0.1, alpha = 0.2) +
          geom_smooth() +
          theme_classic())
  dev.off()
  
}
# }

