library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_MAC1.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/M_MDSC/MDSC_MAC3_MAC1/'
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/M.rds')
cnt <- cnt[, pseudotime]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
cnt <- as.matrix(cnt)
cnt <- cnt[rowMeans(cnt>0.1)>0.01,] ## filter genes
pseudotime = pseudotime[pseudotime %in% cellanno[,1]]
psn <- seq(1,length(pseudotime))
names(psn) <- pseudotime

final <- readRDS('final.rds') 
stat <- data.frame(fdr = final$fdr, foldchange = final$foldchange, stringsAsFactors = F)
stat <- stat[rownames(stat) %in% rownames(cnt), ]
stat <- stat[order(stat$fdr, -stat$foldchange), ]
stat <- stat[stat$fdr<0.05, ]
saveRDS(stat, 'stat.rds')
######### plot
p <- plotGene(testptObj = final,
         Gene = rownames(stat)[1:100],
         Mat = cnt,
         Pseudotime = psn, 
         Design = design,
         Cellanno = cellanno,
         Alpha = 0.2,
         Size = 0.1, FreeScale = T) 
ggsave('top100_gene.png', p, width=22, height=20, dpi = 100)
