library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
# pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_NEU1.rds')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_NEU1.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/MDSC_MAC3_NEU1/'
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/count.rds')
# meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt[, pseudotime]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])
rownames(design) <- as.character(mdsc[,2])
    
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
pseudotime = pseudotime[pseudotime %in% cellanno[,1]]
cnt <- as.matrix(cnt)
rc <- colSums(cnt)
rc <- rc/median(rc)
cnt <- t(log2(t(cnt)/rc + 1))
cnt <- cnt[rowMeans(cnt>1)>0.2,] ## filter genes
final <- readRDS('final.rds') 
res <- final[['res']][rownames(cnt),]
res <- res[order(res[,2]), ]

######### plot
plotGene(Gene = rownames(res)[10:18],
         Mat = cnt,
         Order = data.frame(Cell = pseudotime, Pseudotime = seq(1, length(pseudotime))),
         Design = design,
         Cellanno = cellanno,
         Stat = res,
         Alpha = 0.2,
         Size = 0.1, FreeScale = T)
plotGene(Gene = 'SEPT6',
         Mat = cnt,
         Order = data.frame(Cell = pseudotime, Pseudotime = seq(1, length(pseudotime))),
         Design = design,
         Cellanno = cellanno,
         Stat = res,
         Alpha = 0.2,
         Size = 0.1, FreeScale = F, BySample = TRUE, PlotPoints = TRUE)
