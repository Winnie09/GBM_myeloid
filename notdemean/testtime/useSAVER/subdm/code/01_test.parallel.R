# ------------
# prepare data
# ------------
trajectory <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
# trajectory = 'EMDSC_MAC1'
print(trajectory)
library(parallel)
library(splines)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
library(here)
rdir <- here('notdemean', 'testtime','useSAVER', 'subdm', 'result', trajectory, 'res')
datadir <- here('demean', 'testtime', 'shareEMDSC', 'subdm', 'result', trajectory, 'data') ## data is the same

dir.create(rdir, recursive = TRUE)
pseudotime <- readRDS(paste0(datadir, '/pseudotime.rds'))
expr <- readRDS(paste0(datadir, '/saver.rds'))
cellanno <- readRDS(paste0(datadir, '/cellanno.rds'))
expr <- expr[rowMeans(expr > 0.1) > 0.01, ]

design <- matrix(1, nrow = length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'
# -----
# test
# -----
system.time({
  Res <- testpt(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1.5, verbose=F, ncores=4, type='Time', fit.resolution = 1000, test.pattern = 'overall', demean=FALSE)
  saveRDS(Res, paste0(rdir, '/ptest_res.rds'))
})


