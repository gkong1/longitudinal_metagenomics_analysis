# Package loading 
library(RVAideMemoire) # Functions needed: MVA.synt, MVA.plot, MVA.cmv, MVA.test

## Data loading, filtered and CLR-transformed data (samples as rows, features as columns)
data.comp <- read.delim("tp12_clr_ko_pathways_data.tsv", header = T, row.names = 1)
data.comp <- data.comp[order(rownames(data.comp)),]

metadata <- read.delim(file = "sample_metadata.tsv", header = T, row.names = 1)
metadata$Genotype <- factor(metadata$Genotype, levels = c("WT","HD"))

## Remove outlier
#index <- which(rownames(metadata) == "R5356_S8")
#metadata <- metadata[-index,]
metadata <- metadata[order(rownames(metadata)),]


## Check sample names matched
rownames(data.scaled) == rownames(metadata)

## Set group factor
groupfac <- dummy(metadata$Genotype)

## Run permutation test based on PLS-DA CER
# This may take several mninutes to run, remove progress=FALSE to see computation progress:
MVAtest.res <- MVA.test(X= data.comp, Y = metadata$Genotype, 
                        cmv = FALSE, ncomp = 8, kout = 18, 
                        scale = TRUE, model = c("PLS-DA"))
