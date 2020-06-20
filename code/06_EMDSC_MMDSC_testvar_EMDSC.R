library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/EMDSC_MMDSC_pseudotime.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testvar/E_MDSC/EMDSC_MMDSC/'
dir.create(rdir, showWarnings = F, recursive = T)
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/M.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt[, pseudotime[,1]]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,3])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
    
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
pseudotime = pseudotime[pseudotime[,1] %in% cellanno[,1], ]
cnt <- as.matrix(cnt)
cnt <- cnt[rowMeans(cnt>0.1)>0.01,] ## filter genes

### algo
design = cbind(1, design)
res <- testpt(expr=cnt,cellanno=cellanno,pseudotime=pseudotime,design=design,ncores=8, permuiter=100) 
saveRDS(res, 'final.rds') 
