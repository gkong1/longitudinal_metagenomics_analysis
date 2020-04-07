# Import plasma metabolomics data
library(mixOmics)
mediannormlodata = read.csv("Normalized_metabolites_table.txt", sep = "\t",
                            row.names = 1, header = T, check.names = F)
mediannormlodata = mediannormlodata[order(rownames(mediannormlodata)),] # Sort samples by samplename
mediannormlodata = mediannormlodata[1:18,] # Subset data to experimental groups

# Import metadata
meta = read.csv("Sampleinfo_12.txt", sep = ",", row.names = 4, header = F)
meta = meta[order(rownames(meta)),] # Sort by samplename
colnames(meta) = c("SeqID","Age","Genotype","Cage","SeqRun")
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD"))

# PCA analysis
data.pca = pca(mediannormlodata, ncomp = 10, logratio = "none")
plot(data.pca) # Scree plot
plotIndiv(data.pca, 
          comp = c(1,2), # the components to plot
          pch = meta$Genotype, 
          ind.names = F, 
          group = meta$Genotype, 
          col.per.group = c("#F8766D","#00BFC4"), # Manually set colors
          legend = TRUE,
          title = 'Plasma metabolome at 12wks of age')

# Tuning for number of metabolites to pick for sparse Partial-Least Square Discrimination test.
set.seed(33)  # for reproducible results for this code
data.tune.splsda = tune.splsda(mediannormlodata, 
                                 Y = meta$Genotype, 
                                 ncomp = 3, 
                                 multilevel = NULL, 
                                 logratio = 'none',
                                 test.keepX = c(seq(5,150, 5)), 
                                 validation = c('loo'), 
                                 dist = 'max.dist', 
                                 progressBar = FALSE)
plot(data.tune.splsda)
## Manually set select.keepX according to splsda tuning plot
select.keepX = c(25, 35)

# sPLS-DA 
data.splsda = splsda(X = mediannormlodata,  Y = meta$Genotype, 
                       ncomp = 2, keepX = select.keepX, 
                       logratio= "none")

# sPLS-DA indiv plots
plotIndiv(data.splsda, 
          ind.names = F, 
          col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          X.label= "Component 1",
          Y.label = "Component 2",
          title = 'Plasma metabolome (tp12) sPLS-DA')

## sPLS-DA performance
set.seed(34)  # for reproducible results for this code
data.perf.splsda = perf(data.splsda, validation = 'loo', 
                          progressBar = FALSE, dist = 'max.dist')
data.perf.splsda$error.rate # Shows error rate of sPLS-DA classification model
plot(data.perf.splsda)

## Extract the contribution value and name of the selected metabolite
varcontrib.metab = selectVar(data.splsda, comp = 1)$value
var.comp1.metab = selectVar(data.splsda, comp = 1)$name

## stability of metabolites selected on comp 1
varcontrib.metab$stability = data.perf.splsda$features$stable[[1]][var.comp1.metab]

## Loading plot to visualize importance of each selected metabolite to the classification
plotLoadings(data.splsda, comp = 1, method = 'mean', contrib = 'max',
             title = "Contribution of metabolites on comp 1",
             size.title = 1, ndisplay = 20, size.name = 0.65, size.legend = 0.6)

#--------------------------------------------------------------
# Setting up dataframe for ggplot to plot individual barplots
mediannormlodata1 <- mediannormlodata
mediannormlodata1$Genotype <- factor(meta$Genotype, levels = c("WT","HD"))
mediannormlodata1$Cage <- factor(meta$Cage)
  
## Individual barplot
atp <- ggplot(mediannormlodata1, aes(x=Genotype, y=Adenosine.triphosphate, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("ATP") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

urocanic <- ggplot(mediannormlodata1, aes(x=Genotype, y=Urocanic.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Urocanic acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

carnosine <- ggplot(mediannormlodata1, aes(x=Genotype, y=Carnosine, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Carnosine") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

threonic <- ggplot(mediannormlodata1, aes(x=Genotype, y=Threonic.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Threonic acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

homocitrulline <- ggplot(mediannormlodata1, aes(x=Genotype, y=Homocitrulline, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Homocitrulline") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

orotic <- ggplot(mediannormlodata1, aes(x=Genotype, y=Orotic.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Orotic acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

adp <- ggplot(mediannormlodata1, aes(x=Genotype, y=ADP, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("ADP") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

homocystein <- ggplot(mediannormlodata1, aes(x=Genotype, y=Homocysteine, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Homocysteine") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

paba <- ggplot(mediannormlodata1, aes(x=Genotype, y=p.Aminobenzoic.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("p.Aminobenzoic.acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

trigonelline <- ggplot(mediannormlodata1, aes(x=Genotype, y=Trigonelline, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Trigonelline") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

pipe.acid <- ggplot(mediannormlodata1, aes(x=Genotype, y=Pipecolic.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Pipecolic acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

isobutyrylglycine <- ggplot(mediannormlodata1, aes(x=Genotype, y=Isobutyrylglycine, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Isobutyrylglycine") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

hippurate <- ggplot(mediannormlodata1, aes(x=Genotype, y=Hippuric.acid, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("Hippuric acid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

methylhistidine <- ggplot(mediannormlodata1, aes(x=Genotype, y=Methylhistidine3, fill=Genotype)) + 
  geom_boxplot(aes(fill=Genotype)) + ggtitle("3-Methylhistidine") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  scale_y_continuous(name = "lg") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c('#388ECC','#F68B33'))



