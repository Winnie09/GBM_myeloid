phate <-
  read.csv(
    '/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/Christina_share/m_subsetcelltypes/phate_sub.csv',
    as.is = T
  )

dm <-
  read.csv(
    '/home-4/whou10@jhu.edu/work-zfs/whou10/GBM/data/Christina_share/m_subsetcelltypes/diffusion.csv',
    as.is = T
  )
str(phate)
str(dm)
all.equal(phate[, 1], dm[, 1]) ## exactly the same set of cells

