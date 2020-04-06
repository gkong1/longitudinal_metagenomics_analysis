# Reading in SCFA data
scfa <- read.csv("matrix with FCS_conc.csv", header = T, row.names = 1)

# Statistics with One-way ANOVA
summary(aov(Acetate..uM. ~ Group, data = scfa))
summary(aov(Propionate..uM. ~ Group, data = scfa))
summary(aov(Butyrate..uM. ~ Group, data = scfa))

# Barplots 
library(ggplot2)
acet <- ggplot(scfa, aes(x=Group, y=Acetate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group)) + ggtitle("Acetate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33'))

prop <- ggplot(scfa, aes(x=Group, y=Propionate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group)) + ggtitle("Propionate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33'))

butyr <- ggplot(scfa, aes(x=Group, y=Butyrate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group)) + ggtitle("Butyrate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33'))

gridExtra::grid.arrange(acet, prop, butyr, nrow = 1, ncol = 3)
