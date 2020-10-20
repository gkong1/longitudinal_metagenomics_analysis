# Package loading
library(mixOmics)

# Data loading
counts <- read.delim("alltp_genes_counts.tsv", header = T, row.names = 1)
metadata <- read.delim("sample_metadata.tsv", header = T, row.names = 1)

# Check sample names match
rownames(counts) == rownames(metadata)

# Remove outlier
index <- which(rownames(metadata) == "R5356_S8")
metadata <- metadata[-index,]
counts <- counts[-index,]

# Add pseudocount
counts = counts + 0.1

# Filter low counts and transform to relative abundance
## Function for filtering low counts
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.1 # cutoff chosen
){
  keep.otu = which(rowSums(data)*100/(sum(rowSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

## Actual filter step
result.filter = low.count.removal(counts, percent=0.1)
length(result.filter$keep.otu) 
data.filter = result.filter$data.filter
data.TSS = t(apply(data.filter, 1, TSS.divide))

## Subset to other ages
index = which(meta$Age == "4")
data.TSS.x <- data.TSS[index,]
meta.x = metadata[index,]

rownames(meta.x) == rownames(data.TSS.x)

# PCA analysis
data.pca = pca(data.TSS.x, ncomp = 10, logratio = "CLR")

## Scree plot
plot(data.pca)

## Sample plot
plotIndiv(data.pca, 
          comp = c(1,2), # the components to plot
          ind.names = meta.x$Cage, 
          group = meta.x$Genotype, 
          legend = TRUE) 

# Tuning for sPLS-DA
set.seed(33)
data.tune.splsda = tune.splsda(data.TSS.x, Y = meta.x$Genotype, logratio = 'CLR', ncomp = 4,
                               test.keepX = c(seq(5, 150, 2)), validation = 'loo', 
                               dist = 'centroids.dist', progressBar = F)

plot(data.tune.splsda)
select.keepX = data.tune.splsda$choice.keepX[1:2]
#select.keepX = c(20,25)

# sPLS-DA
data.splsda.x = splsda(data.TSS.x, meta.x$Genotype, ncomp = 2,
                       keepX = select.keepX, logratio = 'CLR')

## Sample plot
plotIndiv(data.splsda.x, col.per.group = color.mixo(1:2),
          comp = c(1,2), pch = meta.x$Genotype, 
          ind.names = F,
          ellipse = T, legend = T)

## Loading plot
plotLoadings(data.splsda.x, comp = 1, method = 'median', contrib = 'max',
             size.title = 1, ndisplay = 7, size.name = 0.7, size.legend = 0.6,
             title = "Contribution on Comp 1")

# Evaluate sPLS-DA performance
set.seed(34)  # for reproducible results for this code
data.perf.splsda = perf(data.splsda.x, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')
data.perf.splsda$error.rate
plot(data.perf.splsda)

# View variables selected by sPLS-DA
head(selectVar(data.splsda.x, comp = 1)$value) 
varcontrib = selectVar(data.splsda, comp = 1)$value
var.comp1 = selectVar(data.splsda, comp = 1)$name
varcontrib$stability = data.perf.splsda$features$stable[[1]][var.comp1]

## Export CLR-transformed data
data.clr = logratio.transfo(data.TSS, logratio = "CLR", offset = 0)
data.good <- as.data.frame(matrix(ncol = ncol(data.TSS), 
                                  nrow = nrow(data.TSS)))
rownames(data.good) <- rownames(data.TSS)
colnames(data.good) <- colnames(data.TSS)
for( i in c(1:nrow(data.TSS))){
  for( j in c(1:ncol(data.TSS))){
    data.good[i,j] <- data.clr[i,j]
  }
}

write.table(data.good, file = "alltp_clr_data.tsv", sep = "\t", quote = F)