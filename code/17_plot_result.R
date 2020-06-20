library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/EMDSC_MMDSC_PMNMDSC_pseudotime.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/EMDSC_MMDSC_PMNMDSC/'
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/M.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt[, pseudotime[,1]]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])  ## 3 is E-MDSC, 8 is M-MDSC
rownames(design) <- as.character(mdsc[,2])
    
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
pseudotime = pseudotime[pseudotime[,1] %in% cellanno[,1],]
cnt <- as.matrix(cnt)
cnt <- cnt[rowMeans(cnt>0.1)>0.01,] ## filter genes
psn <- seq(1, length(pseudotime))
names(psn) <- pseudotime
design = cbind(1, design)
## read result
res <- readRDS(paste0(rdir,'final.rds'))
id = intersect(rownames(cnt), names(res$fdr))
stat <- data.frame(fdr = res$fdr[id], foldchange = res$foldchange[id], stringsAsFactors = FALSE)
stat <- stat[order(stat$fdr, -(stat$foldchange)), ]
stat <- stat[stat$fdr <0.05, ]
write.csv(stat, 'fdr_foldchange.csv')
g = rownames(stat)[1:100]
p <- plotGene(testptObj = res, Gene=g, Mat=cnt[g, ], Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')
ggsave('top100_gene.png', p, width=22, height=20, dpi = 100)

