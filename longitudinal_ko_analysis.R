setwd("~/phd/shotgun/mg-rast/")

## read in KEGG data for all timepoints
kegg_counts <- read.csv("~/phd/shotgun/mg-rast/alltp_mgrast_KO_l3_counts.tsv", sep = "\t", header = T,row.names = 1)
kegg_counts[is.na(kegg_counts)]<- 0
dim(kegg_counts)

kegg_names <- rownames(kegg_counts)
kegg_names1 <- substr(kegg_names, 1, 5)
kegg_names2 <- paste("KO",kegg_names1, sep = "")

rownames(kegg_counts) = kegg_names2
kegg_counts = kegg_counts[,order(colnames(kegg_counts))]

## Reading in metadata
meta <- read.csv("~/phd/shotgun/Sampleinfo_all.csv", sep = ",", header = T, row.names = 1)
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD"))
meta = meta[order(rownames(meta)),]

## Checking whether sample names match
colnames(kegg_counts) == rownames(meta)
colnames(kegg_counts) = rownames(meta)

## Removing outlier
index <- which(rownames(meta) == "R5356_S8")
meta = meta[-index,]
kegg_counts = kegg_counts[-index,]

#----------------------------------------
# Alpha div analysis of raw KO counts
library(phyloseq)
## Importing data to phyloseq
phy.obj.raw = phyloseq(otu_table(kegg_counts, taxa_are_rows = F), sample_data(meta))
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

result.filter = low.count.removal(kegg_counts, percent=0.001)
data.filter = result.filter$data.filter
length(result.filter$keep.otu) # check the number of variables kept after filtering

#-----------------------------
# Calculating bray curtis distance for all tps
library(vegan)
beta_dist_ko <- vegdist(data.filter, index = "bray")
mds_ko <- metaMDS(beta_dist_ko)
mds_data_ko <- as.data.frame(mds_ko$points)

mds_data_ko$Genotype = factor(meta$Genotype, levels = c("WT",'HD'))
mds_data_ko$Age = factor(meta$Age)

# Plotting 
library(ggplot2)
ggplot(mds_data_ko, aes(x=MDS1, y=MDS2)) + 
  geom_jitter(aes(color = Genotype)) + ggtitle("Bray Curtis distance of KOs across timepoints") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right") +
  scale_color_manual(values = c('#388ECC','#F68B33')) +
  facet_grid(~Age) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2")

vegan::adonis(beta_dist_ko ~ Genotype*Age, data = meta)

## Subset data to tp 12
index <- which(mds_data_gene$Age == "12")
mds_data_gene_12 <- mds_data_gene[index,]
meta1_12 = meta[which(meta$Age == "12"),]
vegan::adonis(mds_data_gene_12 ~ Genotype, data = meta1_12)

#-----------------------------------------------------------------------------------
# splsda analysis for tp12 using Mixomics
data.filter = data.filter + 0.00001 # add pseudocount 
## Change to relative abundance
TSS.divide = function(x){
  (x/sum(x))*100
}
data.TSS = t(apply(data.filter, 1, TSS.divide)) # function is applied to each row (i.e. each sample)

## Subset data to tp12 
index = which(meta1$Age == "12")
data.TSS1.12 <- data.TSS1[index,]

## Analysis in Mixomics
library(mixOmics)
### Define number of variables to keep
select.keepX = c(7, 20) # originally c(50,35)

### Actual sPLS-DA
data.splsda = splsda(X = data.TSS1.12,  Y = meta1.12$Genotype, 
                     ncomp = 2, keepX = select.keepX, 
                     logratio= "CLR")

### Sample plot to see how well the samples are distinguished
plotIndiv(data.splsda, 
          ind.names = F, 
          col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta1.12$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          title = 'KEGG (tp12) sPLS-DA')

## Checking performance of sPLS-DA model
set.seed(34)  # for reproducible results for this code
data.perf.splsda = perf(data.splsda, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')
data.perf.splsda$error.rate # Classfication error rate
plot(data.perf.splsda)

# Loading plot to visualize the importance of each gene in their contribution
plotLoadings(data.splsda, comp = 1, method = 'mean', contrib = 'max',
             size.title = 1, ndisplay = 9, size.name = 0.65, size.legend = 0.6)
head(selectVar(data.splsda, comp = 1)$value) 
varcontrib = selectVar(data.splsda, comp = 1)$value # The contribution value by each variable selected
var.comp1 = selectVar(data.splsda, comp = 1)$name # The name of the variable selected

# stability of OTUs selected on comp 1
varcontrib$stability = data.perf.splsda$features$stable[[1]][var.comp1]

#--------------------------------------------------------------
# Individual barplots to examine KO's relative abundances
data.TSS2 <- as.data.frame(data.TSS1)
data.TSS2$Genotype <- factor(meta1$Genotype, levels = c("WT","HD"))
data.TSS2$Age <- factor(meta1$Age)

# Subset data to tp12
data.TSS.12 <- data.TSS2[data.TSS2$Age == "12",]

# Indiv barplot for tp12
## Butanoate metab
buta.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO00650, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Butanoate metabolism") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

## Benzoate degradation
benzo.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO00362 , fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Benzoate degradation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

## Galactose metabolism
galac.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO00052, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Galactose metabolism") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

## Sulfur metabolism
sulf.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO00920, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Sulfur metabolism") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

## Lysine degradation
lysine.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO00310, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Lysine degradation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

## Glutathione metabolism
gluta.12 <- ggplot(data.TSS.12, aes(x=Genotype, y=KO05200, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Glutathione metabolism") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "Relative abundance (%)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

