### This code is for the analysis of bacterial Phylum levels at all timepoints for the longitudinal shotgun sequencing data ###

# Importing Phylum level raw data
count = read.csv("taxa_species_raw_counts.txt", sep = '\t', header = T, row.names = 1)

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

# Change all NA to 0
count[is.na(count)] <- 0

# Removing the rare bacterial species based on raw counts
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
data.filter = low.count.removal(count, percent=0.001)
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
data.TSS = data.TSS[,-index]
index <- which(colnames(data.TSS) == "cannot be assigned to a species ")
data.TSS = data.TSS[,-index] # remove the annoynamous ones

# Analysis in mixOmics
library(mixOmics)

## PCA analysis
data.pca = pca(data.TSS, ncomp = 10, logratio = 'CLR', scale = T, center = T)
plot(data.pca) # Scree plot
plotIndiv(data.pca, 
          comp = c(1,2), 
          group = meta$Genotype, 
          col.per.group = color.mixo(1:2), 
          pch = meta$Cage,
          ind.names = F,
          legend = T,
          title = 'Shotgun species lvl, alltp, PCA comp 1-2')

#----------------------------------------------------------------------------------
# Subset data to tp12
index <- which(meta$Age == "12")
meta.12 <- meta[index,]
data.TSS.12 <- data.TSS[index,]
count.12 <- count[,index]

# Tuning for sPLS-DA at tp12
set.seed(33)
data.tune.splsda = tune.splsda(data.TSS.12, Y = meta.12$Genotype, logratio = 'CLR', 
                               ncomp = 4, test.keepX = c(seq(5, 150, 5)), 
                               validation = 'loo', dist = 'centroids.dist', 
                               progressBar = F)

plot(data.tune.splsda)
select.keepX = c(50,25) # Manually select number of species to keep based on plot(data.tune.splsda)

# sPLS-DA
data.splsda = splsda(data.TSS.12, meta.12$Genotype, ncomp = 2,
                     keepX = select.keepX, logratio = 'CLR')
plotIndiv(data.splsda, 
          ind.names = F, 
          col.per.group = color.mixo(1:2),
          comp = c(1,2), 
          pch = 16, 
          ellipse = T, 
          legend = T, 
          title = 'sPLS-DA shotgun-L7 tp12, comp 1-2')

# Evaluating sPLSDA performance
set.seed(33)
data.perf.splsda = perf(data.splsda, validation = "loo", 
                        progressBar = F, dist = "all")
data.perf.splsda$error.rate #view error rate

# Extract contribution value and name of species selected by sPLS-DA
selectVar(data.splsda, comp = 1)$value
var.comp1 = selectVar(data.splsda, comp = 1)$name

# View stability of species selected by sPLS-DA
data.perf.splsda$features$stable[[1]][var.comp1]

# Plot loadings to visualize contribution of each species in the sPLS-DA model
plotLoadings(data.splsda, comp = 1, method = 'median',
             contrib = 'max', size.title = 1, ndisplay = 10,
             size.name = 0.7, size.legend = 0.7)

#------------------------------------------------------------
# Alpha and beta diversity analysis using Phyloseq
library(phyloseq)

## Alpha div analysis
### Create phyloseq object using raw_counts
phyobj <- phyloseq(otu_table(count.12, taxa_are_rows = TRUE), sample_data(meta.12))

### Rarefy to even depth
set.seed(47)
phy.rare <- rarefy_even_depth(phyobj)

### Test various alpha-div measures
obs = estimate_richness(phy.rare, measures = "Observed")
invs = estimate_richness(phy.rare, measures = "InvSimpson")
shan = estimate_richness(phy.rare, measures = "Shannon")
fisher = estimate_richness(phy.rare, measures = "Fisher")

### Combine to one big dataframe
alpha = as.data.frame(cbind(obs, invs, shan, fisher))
alpha$Genotype = meta.12$Genotype

### Test for significance
kruskal.test(Observed ~ Genotype, data = alpha) 
kruskal.test(Shannon ~ Genotype, data = alpha) 
kruskal.test(Fisher ~ Genotype, data = alpha)
kruskal.test(InvSimpson ~ Genotype, data = alpha)

## Beta-div analysis
### Create phyloseq object
phy.ra = phyloseq(otu_table(data.TSS.12, taxa_are_rows = FALSE), sample_data(meta.12))

### Calculate and ordinate based on Bray Curtis dissimilarity distance
distBC = distance(phy.ra, method = "bray")
ordBC = ordinate(phy.ra, method = "PCoA", distance = distBC)
pcoa_bc = plot_ordination(phy.ra,ordBC, color = "Genotype", 
                          shape = "Genotype") + 
                          geom_point(size = 2, alpha = 0.8) + 
                          ggtitle("PCoA: Bray-Curtis tp12, species lvl") + 
                          theme(plot.title = element_text(hjust = 0.5))

pcoa_bc # Visualize ordination

### Statistical testing with adonis PERMANOVA
vegan::adonis(distBC ~ Genotype, data = meta.12) 


