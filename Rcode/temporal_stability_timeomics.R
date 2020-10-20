#ibrary(shiny)
#runApp('~/Desktop/timeOmics_code/TimeOmics/TimeOmics/')

library(mixOmics)
library(lmtest)
library(dplyr) # for mutate()
library(purrr) # for map()
library(stringr) # for str_split()
library(tseries)
library(tibble) # for rownames_to_columns()
library(tidyr) # for gather()
library(tseries)
library(lmms) 
library(ggplot2)
library(gdata)
library(nlme)
library(reshape2)

# Reading in filtered, CLR-transformed count data
data = read.csv("alltp_clr_counts.txt", sep = "\t", header = T, row.names = 1)

# Reading in metadata
meta = read.csv("sample_metadata.tsv", sep = ',',row.names = 7, header = T)
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD"))

# Check sample names match
rownames(data.good) = rownames(meta)

time_lmms <- rownames(data.good) %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric
sample_id <- rownames(data.good)

# subset WT data
index.wt <- which(meta$Genotype == "WT")
data.good.wt <- data.good[index.wt,]
time_lmms.wt <- time_lmms[index.wt]
sample_id.wt <- rownames(data.good.wt)

# subset HD data
index.hd <- which(meta$Genotype == "HD")
data.good.hd <- data.good[index.hd,]
time_lmms.hd <- time_lmms[index.hd]
sample_id.hd <- rownames(data.good.hd)

###-----------------------------------------------------------------------------------------------
# lmmSpline calculation
### For WT data
data.good.x <- data.good.wt
time_lmms.x <- time_lmms.wt
sample_id.x <- sample_id.wt

spline.cubicpspline.wt = lmmSpline(data = data.good.x, time = time_lmms.x, 
                                sampleID = sample_id.x,
                                basis = 'cubic p-spline', keepModels = T , numCores = 4)
spline.cubicpspline.wt@modelsUsed 


### to get the models used, 0 means straight lines,
### 1-> 3 means less -> more complex curves
# filter splines based on homoskedasticity
## wrapper.filter.splines function in filter.R
load("~/phd/timeOmics_code/filter.R")
## get_MSE function from filter.R

# Filter splines on WT data
filter.cubicpsplines.res.wt <- wrapper.filter.splines(data.good.x, spline.cubicpspline.x, stationnarity.test = F,
                                                      homoskedasticity = T, MSE.filter = F)
index.filter.wt <- (rownames(spline.cubicpspline.x@predSpline) %in% filter.cubicpsplines.res.x$to_keep) %>% which()

# Subset spline data with the ones that passed threshold
spline.data.wt <- spline.cubicpspline.wt@predSpline[index.filter.wt,] %>% t %>% as.data.frame()

model.evo.wt <- spline.data.wt %>% rownames_to_column("time") %>% gather(Features,value, -time) %>% mutate(time=as.numeric(time)) %>% 
    ggplot(aes(x=time, y = value, col = Features)) +
    geom_line() + theme_bw() + theme(legend.position = "none") +
    ggtitle("Modelled OTU evolution on WT data")
model.evo.wt

# Use mixomics:PCA to cluster timecourse profiles
library(mixOmics)

# Perform clustering
## pca.get_cluster() and pca.plot() in tune.spca.R
pca.res.wt <- pca(spline.data.wt, ncomp = 2, scale = T, center = T)
pca.wt.s <- pca.plot(pca.res.wt, title = "PCA, WT data, scale = T") # to plot the clusters
pca.get_cluster(pca.res.wt) %>% pull(cluster) %>% table
plotIndiv(pca.res.wt)
plotVar(pca.res.wt)

# Wrapper.silhouette.pca() in tune.spca.R
# To get silhouette coefficient for this clustering
wrapper.silhouette.pca(spline.data.wt, ncomp = 2, scale = T, center = T)


# Timecourse sparse PCA to clean clusters
source("silhouette.R")
source("tune.spca.R")

keepX = list(seq(11,29,3), seq(9,15,1)) 

# tuning the silhouette criterion to tune clustering and
# to get optimal number of molecule (OTU) per cluster
res.tune.spca.wt <- tune.spca(X = spline.data.wt, ncomp = 2, keepX = keepX)
tune.spca.choice.keepX(res.tune.spca.wt, draw = T)

# Actual sPCA
keepX = c(11,10)
spca.res.wt <- spca(spline.data.wt, ncomp = 2, keepX = keepX) # how to select number of params to keep
pca.get_cluster(spca.res.wt) %>% pull(cluster) %>% table
wrapper.silhouette.spca(spline.data.wt, keepX = c(11,10),
                        ncomp = 2, scale = T, center=T,
                        plot.t=T)


spca.plot(spca.res.wt, "sparse PCA WT data")

selectVar(spca.res.wt, comp = 1) # to obtain the name of variables in spca

spca_clus_wt <- pca.get_cluster(spca.res.wt)
spca_clus_wt <- spca_clus_wt[spca_clus_wt$cluster != "0",]
spca_clus_wt <- spca_clus_wt[order(spca_clus_wt$cluster),]

write.table(spca_clus_wt, file="spca_wt_cluster.txt", sep = "\t", quote=F,row.names = F)
