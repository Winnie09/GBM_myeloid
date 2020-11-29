rm(list = ls())
library(ggplot2)
library(reshape2)
library(RColorBrewer)
suppressMessages(library(igraph))
library(parallel)
library(splines)
library(viridis)
library(here)

source(
  '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R'
)
trajectory = 'EMDSC_MAC1'
# for (trajectory in c('EMDSC_MAC1','EMDSC_MAC1_MAC2','EMDSC_MMDSC','EMDSC_MMDSC_PMNMDSC','EMDSC_MMDSC_MAC1')) {
print(trajectory)
rdir <-
  here('demean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'result',
       trajectory,
       'res')
datadir <-
  here('demean',
       'testtime',
       'shareEMDSC',
       'subdm',
       'result',
       trajectory,
       'data')
plotdir <-
  here('demean', 'testtime', 'shareEMDSC', 'subdm', 'plot', trajectory)

Res <- readRDS(paste0(rdir, '/ptest_res.rds'))
res = data.frame(
  fdr = Res$fdr,
  foldchange = Res$foldchange,
  pvalue = Res$pvalue
)
res = res[res[, 1] < 0.05,]
res = res[order(res[, 1], res[, 2]),]

pt = Res$pseudotime
pred <- Res$predict.values[rownames(res), names(pt)]
str(pred)

# cluster genes
cm <- rowMeans(pred)
csd <- apply(pred, 1, sd)
pred <- (pred - cm) / csd
set.seed(12345)
clu <- kmeans(pred, 5)$cluster
clu <- sort(clu)
table(clu)
res$clu = clu[rownames(res)]

# annotate rows and columns
colann <-
  data.frame(sample = cellanno[match(cellanno[, 1], colnames(pred)), 2],
             celltype = active.ident[colnames(pred)],
             stringsAsFactors = F)
rownames(colann) = colnames(pred)
rowann = data.frame(cluster = as.character(clu), stringsAsFactors = F)
rownames(rowann) = names(clu)

# define colors
pd <- Res$expr.demean[names(clu), names(pt)]
v = loess(pd[1, ] ~ seq(1, ncol(pd)))$fitted
plot(pd[1, ] ~ seq(1, ncol(pd)), pch = 20)
lines(v ~ seq(1, ncol(pd)), col = 'red')
plot(v ~ seq(1, ncol(pd)), col = 'red')

pd <- Res$predict.values[names(clu), names(pt)]
v = loess(pd[1, ] ~ seq(1, ncol(pd)))$fitted
plot(pd[1, ] ~ seq(1, ncol(pd)), pch = 20)
lines(v ~ seq(1, ncol(pd)), col = 'red')
plot(v ~ seq(1, ncol(pd)), col = 'red')

fit.method <- sapply(names(clu), function(i) {
  get_population_fit(testobj = Res,
                     variable = NA,
                     gene = i)
})

pd = Res$predict.values[names(clu),]
fit.pred <- sapply(1:nrow(pd), function(i) {
  loess(pd[i, ] ~ seq(1, ncol(pd)))$fitted
})

pd = Res$expr.demean[names(clu),]
fit.demean <- sapply(1:nrow(pd), function(i) {
  loess(pd[i, ] ~ seq(1, ncol(pd)))$fitted
})

c <-
  sapply(1:ncol(fit.pred), function(i)
    cor(fit.method[, i], fit.pred[, i]))
c <-
  sapply(1:ncol(fit.pred), function(i)
    cor(fit.demean[, i], fit.pred[, i]))


for (i in 1:ncol(fit.pred)) {
  pdf(
    paste0(plotdir, '/', names(clu)[i], '_pattern.pdf'),
    width = 10,
    height = 3
  )
  par(mfrow = c(1, 3))
  plot(
    fit.pred[, i] ~ seq(1, nrow(fit.pred)),
    main = names(clu)[i],
    xlab = 'Pseudotime',
    ylab = 'Fit on Predict.values',
    pch = 20
  )
  plot(
    fit.demean[, i] ~ seq(1, nrow(fit.demean)),
    main = names(clu)[i],
    xlab = 'Pseudotime',
    ylab = 'Fit on Demean Expression',
    pch = 20
  )
  plot(
    fit.method[, i] ~ seq(1, nrow(fit.method)),
    main = names(clu)[i],
    xlab = 'Pseudotime',
    ylab = 'Phi * X * Beta',
    pch = 20
  )
  dev.off()
}

