library(Matrix)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
# pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC3_NEU1.rds')
pseudotime <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/data/order/MDSC_MAC1.rds')
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/result/MDSC_MAC1/'
dir.create(rdir, showWarnings = F, recursive = T)
setwd(rdir)

cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/count.rds')
meta <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/GBM/singleObject/M/meta.rds')
cnt <- cnt[, pseudotime]
cellanno <- data.frame(cell = colnames(cnt), sample = sapply(colnames(cnt), function(i) sub('_.*','',sub('.*-','',i))), stringsAsFactors = FALSE)
mdsc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/GBM/meta/mdsc_proportions.csv', header = T)
design <- data.frame(MdscProp = mdsc[,8])
rownames(design) <- as.character(mdsc[,2])
    
cellanno <- cellanno[cellanno[,2] %in% rownames(design),]
cnt <- cnt[, cellanno[,1]]
pseudotime = pseudotime[pseudotime %in% cellanno[,1]]
cnt <- as.matrix(cnt)
### 
design = cbind(1,design)
library(gtools)
# per <- permutations(nrow(design),nrow(design))
per <- t(sapply(1:100, function(i){
  set.seed(i)
  sample(rownames(design))
}))

id <- which(!duplicated(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))))
per <- per[id,]
oriid <- which(apply(per,1,function(i) paste0(as.vector(design[i,]),collapse = '_'))==paste0(as.vector(design),collapse = '_'))
library(parallel)
# perll <- mclapply(1:2,function(i) {
perll <- mclapply(1:nrow(per),function(i) {
  perdesign <- design[per[i,],,drop=F]
  row.names(perdesign) <- row.names(design)
  diffpt(expr=cnt,design=perdesign,pseudotime=pseudotime,num.knot = 3,cellanno = cellanno, verbose = T)$logL
}, mc.cores =4)
perll <- do.call(cbind,perll)
pval <- sapply(1:nrow(perll), function(i) pnorm(perll[i,oriid],mean(perll[i,-oriid]),sd(perll[i,-oriid]),lower.tail = F))
fdr <- p.adjust(pval,method='fdr')
names(fdr) <- row.names(perll)

saveRDS(fdr, 'fdr.rds')  
res <- data.frame(P.Value = pval, adj.P.Val = fdr, stringsAsFactors = F)
rownames(res) <- names(fdr)
res <- res[order(res[,2]),]
sensfdr <- SensFdr(Order = rownames(res), TruePositive = selgene, statistics=res)
final <- list()
final[['res']] <- res
final[['sensfdr']] <- c(method, AreaUnderSensFdr(sensfdr))
final[['perll']] <- perll
saveRDS(final, 'final.rds')  

