# ------------
# prepare data
# ------------
trajectory <- as.character(commandArgs(trailingOnly = TRUE)[[1]][1])
variable <- as.character(commandArgs(trailingOnly = TRUE)[[2]][1])
# trajectory = 'EMDSC_MAC1'
# variable = 'EMDSC'
library(parallel)
library(splines)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/', trajectory, '/', variable, '/')
datadir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/', trajectory, '/data/')
pseudotime <- readRDS(paste0(datadir, 'pseudotime.rds'))
expr <- readRDS(paste0(datadir, 'log2cpm.rds'))
cellanno <- readRDS(paste0(datadir, 'cellanno.rds'))
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv')

rmcell <- cellanno[cellanno[,2]=='GBM057',1]
expr <- expr[, !colnames(expr)%in%rmcell]
pseudotime <- pseudotime[!names(pseudotime) %in% rmcell]
cellanno <- cellanno[!cellanno[,2]%in%'GBM057',]


if (variable == 'EMDSC'){
  design <- data.frame(intercept = 1, emdscProp = mdsc[,3])
} else if (variable == 'MMDSC'){
  design <- data.frame(intercept = 1, mmdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
} else {
  print('variable should be either EMDSC or MMDSC!')
}

rownames(design) <- as.character(mdsc[,2])

# -----
# test
# -----
system.time({
  Res <- ptest(expr = expr, cellanno = cellanno, pseudotime = pseudotime, design=design, permuiter=100, EMmaxiter=100, EMitercutoff=1, verbose=F, ncores=1, type='Variable', fit.resolution = 1000)
  saveRDS(Res, paste0(rdir, 'ptest_res.rds'))
})

