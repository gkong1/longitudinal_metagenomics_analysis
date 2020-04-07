### This code is for the analysis of bacterial Phylum levels at all timepoints for the longitudinal shotgun sequencing data ###

# Importing Phylum level raw data
counts = read.csv("taxa_phylum_raw_counts.txt", sep = "\t", header = T, row.names =1)

# Importing metadata
meta <- read.csv("Sampleinfo_all.csv", sep = "\t", header = F, row.names = 1)
meta = meta[order(rownames(meta)),] # Sort by samplename
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD")) # Set Genotype factor levels
meta$Age = factor(meta$Age, levels = c("4","6","8","10","12")) # Set Age factor levels

rownames(meta) == colnames(counts) # Check whether samples are in same order in both dataframe

# Remove Outlier "R5356_S8"
index <- which(rownames(meta) == "R5356_S8")
meta = meta[-index,]
counts = counts[,-index]

# Removing the rare bacterial Phylum based on raw counts
## Function for low count removal
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.001 # cutoff chosen
){
  keep.otu = which(rowSums(data)*100/(sum(rowSums(data))) > percent)
  data.filter = data[keep.otu,]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}
## Filtering
data.filter = low.count.removal(counts)
data.filter = data.filter$data.filter
data.filter = data.filter + 1 # Add pseudocount

# Convert to relative abundance
## Function
TSS.divide = function(x){
  (x/sum(x))*100
}
## Tranforming to relative abundance
data.TSS = t(apply(data.filter, 2, TSS.divide))

## Removed the unclassified
index <- which(colnames(data.TSS) == "unclassified")
data.TSS <-  data.TSS[,-index]
index <- which(colnames(data.TSS) == "cannot be assigned to a phylum ")
data.TSS <- data.TSS[,-index]

# Analysis in mixOmics
library(mixOmics)

## PCA analysis
data.pca = pca(data.TSS, ncomp = 10, logratio = 'CLR')
plot(data.pca)
plotIndiv(data.pca, 
          comp = c(1,2), 
          pch = meta$Age, 
          ind.names = F, 
          group = meta$Genotype, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title = 'PCA comp 1 - 2')

# Individual barplots of most abundant bacterial Phylum
## Sort by phylum abundance
index = order(colSums(data.TSS),decreasing = T)
data.TSS_1 <- data.TSS[,order(colSums(data.TSS),decreasing = T)]
data.TSS_1 <- data.TSS[,order(colSums(data.TSS),decreasing = T)]
data.TSS_1$Genotype = factor(meta$Genotype, levels = c("WT","HD"))

## Indiv Barplots of top 4 most abundant Phylum
bacteroidetes <- ggplot(data.TSS_1, aes(x=Genotype, y=Bacteroidetes, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Bacteroidetes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  facet_grid(. ~ Age) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

firmicutes <- ggplot(data.TSS_1, aes(x=Genotype, y=Firmicutes, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Firmicutes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  facet_grid(. ~ Age) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

proteobacteria <- ggplot(data.TSS_1, aes(x=Genotype, y=Proteobacteria, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Proteobacteria") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  facet_grid(. ~ Age) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

actinobacteria <- ggplot(data.TSS_1, aes(x=Genotype, y=Actinobacteria, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Actinobacteria") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  facet_grid(. ~ Age) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))
