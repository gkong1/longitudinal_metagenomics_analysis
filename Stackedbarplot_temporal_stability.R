setwd("~/phd/02_Refined_analyses/")
data <- read.csv("Temporal_profile_summary.txt", sep = "\t", header = T)
data$Type = factor(data$Type, levels = c("Species","KO","Genes"))
data$Genotype = factor(data$Genotype, levels = c("WT","HD"))

library(ggplot2)
library(viridis)

# Stacked barplots
jpeg("Temporal profile cluster analysis.jpeg", width = 6.5, height = 5.7, units = "in", res = 300)
ggplot(data, aes(fill= Model, y= Value, x= Genotype)) + 
  geom_bar(position="fill", stat="identity") + facet_grid(~Type) + 
  scale_fill_viridis(discrete = T) +
  ylab("Proportion") +
  ggtitle("Temporal profile cluster analysis") 
dev.off()
