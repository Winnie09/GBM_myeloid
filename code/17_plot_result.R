library(Matrix)
path = 'EMDSC_MMDSC_PMNMDSC'
# path = 'EMDSC_MMDSC_MAC1'
# path = 'EMDSC_MMDSC'
# path = 'EMDSC_MAC1'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pseudotime <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/', path, '_pseudotime.rds'))
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/testtime/M_MDSC/', path, '/')
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



#####  use intersect DEG to doclustering
### get fitted values
fit = get_fitted_values(testObj = res, Gene = rownames(stat), Pseudotime = pseudotime[,2])
fit_scale <- t(apply(fit, 1, scale))

### elbow's method, determine number of clsuters
bttss <- sapply(1:20, function(i){
  set.seed(12345)
  clures = kmeans(fit_scale,  centers = i)
  (clures$totss - sum(clures$withinss))/clures$totss
})
plot(bttss ~ seq(1,20), pch=20, xlab='Number of clusters', ylab='BTSS/TOTSS')

### clustering
set.seed(12345)
clu = kmeans(fit_scale,  centers = 5)$cluster
df = data.frame(stat, Cluster = clu[rownames(stat)])
write.csv(df, 'fdr_foldchange_cluster.csv')

### plot
par(mfrow=c(2,5))
for (i in  seq(1, max(clu))){
  plot(colMeans(fit_scale[names(clu[clu==i]), ,drop=F]), pch = 20, xlab='Pseudotime', ylab='Scaled Expression', main = paste0('clu', i, ':', length(clu[clu==i])))
}

selectgene = c(names(clu[clu==1][1]), names(clu[clu==2][1]), names(clu[clu==3][1]), names(clu[clu==4][1]), names(clu[clu==5][1]))
for (g in selectgene){
  plot(fit[g,], pch = 20, xlab='Pseudotime', ylab='Expression', main = g)  
}


p1 <- plotGene(testptObj = res, Gene=selectgene[1], Mat=cnt, Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time') + theme(legend.position='none')


p2 <- plotGene(testptObj = res, Gene=selectgene[2], Mat=cnt, Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')+ theme(legend.position='none')


p3 <- plotGene(testptObj = res, Gene=selectgene[3], Mat=cnt, Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')+ theme(legend.position='none')


p4 <- plotGene(testptObj = res, Gene=selectgene[4], Mat=cnt, Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')+ theme(legend.position='none')


p5 <- plotGene(testptObj = res, Gene=selectgene[5], Mat=cnt, Pseudotime = pseudotime, Cellanno=cellanno, Design=design,  Alpha=0.5, Size=0.1, PlotPoints = F, FreeScale = TRUE, BySample = F, type = 'Time')+ theme(legend.position='none')
library(gridExtra)
grid.arrange(p1,p2,p3,p4,p5, nrow=2)



