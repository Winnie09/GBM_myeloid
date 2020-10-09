library(destiny)
cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/M.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
active.ident <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/active.ident.rds')
meta = cbind(meta, active.ident = active.ident[rownames(meta)])
meta <- meta[meta$Grade == 'High' & meta$Treatment == 'Untreated', ]
selectcell <- rownames(meta[meta$active.ident %in% c('MDSC', 'MAC1', 'MAC2', 'MAC3', 'NEU1', 'MAC4'), ])
cnt <- cnt[, selectcell]
saveRDS(cnt, '/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/6celltype_SAVER.rds')
cnt <- cnt[rowMeans(cnt>0.1)>0.01, ]

cv = apply(cnt, 1, sd)/rowMeans(cnt)
dm2 <- DiffusionMap(cnt[names(sort(cv, decreasing = T)[1:3e3]),])
saveRDS(dm2, 'myeloid_6ct_3e3cvgene_cv.rds')

