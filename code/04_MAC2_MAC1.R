library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
# pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_NEU1.rds')
# pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC1.rds')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MAC2_MAC1.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/M_MDSC/MAC2_MAC1/'
dir.create(rdir, showWarnings = F, recursive = T)
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/count.rds')
# # ### subset a small test set
# set.seed(12345)
# id1 = sample(rownames(cnt),3)
# cnt <- cnt[id1, ]
##### 
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt[, pseudotime]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
    
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
pseudotime = pseudotime[pseudotime %in% cellanno[,1]]
cnt <- as.matrix(cnt)
rc <- colSums(cnt)
rc <- rc/median(rc)
cnt <- t(log2(t(cnt)/rc + 1))
cnt <- cnt[rowMeans(cnt>0.1)>0.01,] ## filter genes

### algo
psn <- seq(1, length(pseudotime))
names(psn) <- pseudotime
design = cbind(1, design)
res <- testpt(expr=cnt,cellanno=cellanno,pseudotime=psn,design=design,ncores=1, permuiter=100) 
saveRDS(res, 'final.rds') 
