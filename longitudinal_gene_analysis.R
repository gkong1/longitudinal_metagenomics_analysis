## Importing gene_data for all tp
gene_counts <- read.csv("~/phd/shotgun/mg-rast/alltp_mgrast_KO_function_counts_names.txt", sep = "\t", header = T)
gene_names = gene_counts[,1] # used in the case to keep original names
rownames(gene_counts) = make.names(gene_names, unique = T) ## Problem: has put X in front of some of the gene names
gene_counts = gene_counts[,-1] 
gene_counts[is.na(gene_counts)]<- 0
dim(gene_counts)

gene_counts = gene_counts[,order(colnames(gene_counts))] # Sample by number order
gene_counts = t(gene_counts) # To have samples as rows x genes as cols

## Reading in metadata
meta <- read.csv("~/phd/shotgun/Sampleinfo_all.csv", sep = ",", header = T, row.names = 1)
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD"))
meta = meta[order(rownames(meta)),]

rownames(gene_counts) = rownames(meta) # Making sure sample names are the same

## Removing outlier
index <- which(rownames(meta) == "R5356_S8")
meta = meta[-index,]
gene_counts = gene_counts[-index,]

#-----------------------------
# Alpha div analysis of raw gene counts
library(phyloseq)
## Importing data to phyloseq
phy.obj.raw = phyloseq(otu_table(gene_counts, taxa_are_rows = F), sample_data(meta))
## Rarefying data to min library size
set.seed(223)
phy.obj.rare = rarefy_even_depth(phy.obj.raw, sample.size = min(sample_sums(phy.obj.raw)))
## Alpha div analysis
alpha_div = plot_richness(phy.obj.rare, shape = "Genotype", color = "Genotype", x = "Age",
                          measures = c("Observed", "Shannon", "InvSimpson","Fisher"), 
                          title = "Alpha Diversity")

#-----------------------------------
# Low abundance filtering
## Function
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.1 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

## Perform filtering
result.filter = low.count.removal(gene_counts, percent=0.1)
data.filter = result.filter$data.filter
length(result.filter$keep.otu) # check the number of variables kept after filtering

#-----------------------------
# Calculating bray curtis distance for all tps
library(vegan)
beta_dist_gene <- vegdist(data.filter, index = "bray")
mds_gene <- metaMDS(beta_dist_gene)
mds_data_gene <- as.data.frame(mds_gene$points)

mds_data_gene$Genotype = factor(meta$Genotype, levels = c("WT",'HD'))
mds_data_gene$Age = factor(meta$Age)

# Plotting 
library(ggplot2)
ggplot(mds_data_gene, aes(x=MDS1, y=MDS2)) + 
  geom_jitter(aes(color = Genotype)) + ggtitle("Bray Curtis distance of genes across timepoints") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right") +
  scale_color_manual(values = c('#388ECC','#F68B33')) +
  facet_grid(~Age) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2")

vegan::adonis(beta_dist_gene ~ Genotype*Age, data = meta)

## Subset data to tp 12
index <- which(mds_data_gene$Age == "12")
mds_data_gene_12 <- mds_data_gene[index,]
meta1_12 = meta[which(meta$Age == "12"),]
vegan::adonis(mds_data_gene_12 ~ Genotype, data = meta1_12)

#---------------------------------
# sPLS-DA analysis
## Add pseudocount for mixomics
data.filter = data.filter + 0.000001
## Change to relative abundance
TSS.divide = function(x){
  (x/sum(x))*100
}
data.TSS = t(apply(data.filter, 1, TSS.divide)) # function is applied to each row (i.e. each sample)

## Subset data to tp 12
index = which(meta$Age == "12")
data.TSS.12 <- data.TSS[index,]

select.keepX = c(20, 20) # selecting number of desired variables

## Actual sPLS-DA
data.splsda.12 = splsda(X = data.TSS.12,  Y = meta.12$Genotype, 
                     ncomp = 2, keepX = select.keepX, 
                     logratio= "CLR")
## Sample plot to see how well the samples are distinguished
plotIndiv(data.splsda.12, 
          ind.names = F, 
          col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta.12$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          X.label = "Component 1",
          Y.label = "Component 2",
          title = 'Genes (tp12) sPLS-DA')

## Checking performance of sPLS-DA model
set.seed(34)  # for reproducible results for this code
data.perf.splsda = perf(data.splsda.12, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')
data.perf.splsda$error.rate # Classfication error rate
plot(data.perf.splsda)

# Loading plot to visualize the importance of each gene in their contribution
plotLoadings(data.splsda.12, comp = 1, method = 'median', contrib = 'max',
             size.title = 1, ndisplay = 15, size.name = 0.7, size.legend = 0.6,
             title = "Contribution of genes on Comp 1")

head(selectVar(data.splsda, comp = 1)$value) 
varcontrib = selectVar(data.splsda, comp = 1)$value # The contribution value by each variable selected
var.comp1 = selectVar(data.splsda, comp = 1)$name # The name of the variable selected

# stability of genes selected on comp 1
varcontrib$stability = data.perf.splsda$features$stable[[1]][var.comp1]

#--------------------------------
# Individual barplots to examine gene relative abundances at tp12
data.TSS2 <- as.data.frame(data.TSS.12)
data.TSS2$Genotype <- factor(meta.12$Genotype, levels = c("WT","HD"))
data.TSS2$Age <- factor(meta.12$Age)

gmd <- ggplot(data.TSS2, aes(x=Genotype, y=E4.2.1.47, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("gmd") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

pyrE <- ggplot(data.TSS2, aes(x=Genotype, y=pyrE, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("pyrE") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

pepD <- ggplot(data.TSS2, aes(x=Genotype, y=K01270, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("pepD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

oadB <- ggplot(data.TSS2, aes(x=Genotype, y=E4.1.1.3B, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("oadB") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

ara <- ggplot(data.TSS2, aes(x=Genotype, y=araA, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("araA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

pyrG <- ggplot(data.TSS2, aes(x=Genotype, y=E6.3.4.2, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("pyrG") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))


