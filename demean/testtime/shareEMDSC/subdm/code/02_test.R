# ------------
# prepare data
# ------------
trajectory <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
# trajectory = 'EMDSC_MAC1'
print(trajectory)
library(parallel)
library(splines)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/res/')
datadir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/data/')
dir.create(rdir, recursive = TRUE)
pseudotime <- readRDS(paste0(datadir, 'pseudotime.rds'))
expr <- readRDS(paste0(datadir, 'log2cpm.rds'))
cellanno <- readRDS(paste0(datadir, 'cellanno.rds'))

expr <- expr[rowMeans(expr > 0) > 0.01, ]
design <- matrix(1, nrow = length(unique(cellanno[,2])))
rownames(design) <- unique(cellanno[,2])
colnames(design) <- 'intercept'
## demean
expr.demean <- lapply(unique(cellanno[,2]), function(s){
  tmp <- expr[, cellanno[cellanno[,2] == s, 1]]
  tmp2 <- tmp- rowMeans(tmp)
})
expr.demean <- do.call(cbind, expr.demean)

# -----
# test
# -----
system.time({
  Res <- testpt(expr = expr.demean, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=12, type='Time', fit.resolution = 1000, test.pattern = 'overall')
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})


