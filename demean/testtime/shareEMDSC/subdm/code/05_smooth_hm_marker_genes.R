rm(list=ls())
library(ggplot2)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')

## read data
trajectory = 'EMDSC_MMDSC'
print(trajectory)
rdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/result/', trajectory, '/res/')
plotdir <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/GBM_myeloid/demean/testtime/subdm/plot/', trajectory, '/')

Res <- readRDS(paste0(rdir, 'ptest_res.rds'))
pt = Res$pseudotime

res <- read.csv(paste0(rdir, 'differential_genes_clu.csv'), row.names = 1)
clu = res$clu
names(clu) = rownames(res)
pred <- Res$predict.values[rownames(res), names(pt)]

gl <- c('LYZ', 'EREG', 'TIMP1', 'S100A6', 'CXCL2', 'CXCL3', 'IL1B', 'MARCO', 'CCL20', 'PTGER2', 'CD44', 'BCL2A1', 'AREG', 'S100A9', 'ANXA1', 'CD93', 'HK2', 'SLC2A5', 'ENO2', 'GPNMB', 'HMOX1', 'SLC2A1')

# ## fit using GAM on original expression
# fit <- sapply(1:nrow(Res$expr.ori), function(i){
#   model <- mgcv::gam(Res$expr.ori[i,]~s(pt,k=3))
#   fitted(model)
# })
# p <- STIP(fit, gl)  

## fit using loess on predicted (sample-level) expression
fit <- t(sapply(1:nrow(pred), function(i){
  tmp <- data.frame(expr = pred[i,names(pt)], 
                    pseudotime = pt,
                    stringsAsFactors = F)
  v <- loess(expr~pseudotime, data = tmp, degree = 2)$fitted
}))
rownames(fit) <- rownames(pred)
colnames(fit) <- seq(1, ncol(fit))
saveRDS(fit, paste0(rdir, 'gene_by_cell_smooth_fit.rds'))

## apply STIP to plot
STIP <- function(fit,gl) {
      dn <- dimnames(fit)
      fit <- t(apply(fit,1,scale))
      dimnames(fit) <- dn
      gene <- row.names(fit)
      
      zpdirection <- fit[,1] < fit[,ncol(fit)]
      
      zp <- apply(fit,1,function(sf) {
            names(which(sapply(1:(length(sf)-1),function(i) sf[i]*sf[i+1] < 0)))
      })
      zpnum <- sapply(zp,length)
      inczp <- names(which(zpdirection[zpnum==1]))
      deczp <- names(which(!zpdirection[zpnum==1]))
      multipoint <- names(zpnum)[zpnum > 1]
      m1 <- names(which(fit[multipoint,1] > 0))
      m2 <- names(which(fit[multipoint,1] < 0))
      
      geneorder <- NULL
      
      if (length(deczp) > 0) {
        tmp <- unlist(zp[deczp])
        n <- names(tmp)
        tmp <- as.numeric(tmp)
        names(tmp) <- n
            geneorder <- c(geneorder,names(sort(tmp,decreasing = F)))
      }
      geneorder <- c(geneorder,names(sort(sapply(zp[m2],function(i) as.numeric(i[1])))))
      if (length(inczp) > 0) {
        tmp <- unlist(zp[inczp])
        n <- names(tmp)
        tmp <- as.numeric(tmp)
        names(tmp) <- n
            geneorder <- c(geneorder,names(sort(tmp)))
      }
      geneorder <- c(geneorder,names(sort(sapply(zp[m1],function(i) as.numeric(i[1])))))
      geneorder <- rev(geneorder)
      plotdata <- fit[geneorder,]
      plotdata <- melt(plotdata)
      colnames(plotdata) <- c("Gene","Pseudotime","Expression")
      
      p1 <- ggplot(plotdata,aes(Pseudotime,Gene,fill=Expression)) + geom_tile() + theme_classic() + scale_fill_gradient2(low="blue",high="red",midpoint=0) + theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),plot.margin=margin(2,0,2,0)) + scale_x_continuous(expand=c(0.001,0.001))
      
      yax <- rep("A",length(gene))
      inid <- which(gl %in% gene)
      oid <- which(!gl %in% gene)
      yaxglid <- round(seq(1,length(gene)-2,length.out=length(inid)))
      yax[yaxglid] <- gl
      yax[setdiff(1:length(yax),yaxglid)] <- setdiff(gene,gl)
      p2 <- ggplot() + geom_point(data=data.frame(gene=factor(yax,levels=yax),x=1),aes(x=x,y=gene),col="white") + geom_text(data=data.frame(text=gl,id=gl,x=1),aes(x=x,y=id,label=gl),size=3) + geom_segment(data=data.frame(x=1.5,xend=2,y=gl,yend=yax[match(gl,geneorder)]),aes(x=x,y=y,xend=xend,yend=yend),size=0.1) + theme_classic() + coord_cartesian(xlim=c(0.5,2)) + theme(axis.title = element_text(color="white"),axis.line = element_line(color="white"),axis.text = element_text(color="white"),axis.ticks = element_line(color="white"),plot.margin=margin(2,0,2,2)) + scale_x_continuous(expand=c(0,0))
      gridExtra::grid.arrange(p2,p1,nrow=1,layout_matrix=cbind(1,2,2,2))
}

pdf(paste0(plotdir, 'smooth_hm_mark_genes.pdf'), width = 8, height = 9)  
STIP(fit, gl)  
dev.off()

