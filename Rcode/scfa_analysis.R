
# Reading in SCFA data
scfa <- read.csv("matrix with FCS_conc.csv", header = T, row.names = 1)

# F-statistics with var.test to test for equal variance
var.test(Acetate..uM. ~ Group, data = scfa)
var.test(Propionate..uM. ~ Group, data = scfa)
var.test(Butyrate..uM. ~ Group, data = scfa)

# Statistics with Welch T-Test for unequal variance
t.test(Acetate..uM. ~ Group, data = scfa, var.equal = FALSE)
t.test(Propionate..uM. ~ Group, data = scfa, var.equal = FALSE)
t.test(Butyrate..uM. ~ Group, data = scfa, var.equal = FALSE)

# Statistics with One-way ANOVA
summary(aov(Acetate..uM. ~ Group, data = scfa))
summary(aov(Propionate..uM. ~ Group, data = scfa))
summary(aov(Butyrate..uM. ~ Group, data = scfa))

# Barplots 
library(ggplot2)
acet <- ggplot(scfa, aes(x=Group, y=Acetate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group), outlier.shape = NA) + ggtitle("Acetate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33')) + geom_jitter()

prop <- ggplot(scfa, aes(x=Group, y=Propionate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group), outlier.shape = NA) + ggtitle("Propionate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33')) + geom_jitter()

butyr <- ggplot(scfa, aes(x=Group, y=Butyrate..uM., fill=Group)) + 
  geom_boxplot(aes(fill=Group), outlier.shape = NA) + ggtitle("Butyrate") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(element_blank()) +
  ylab("Concentration (uM)") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c('#388ECC','#F68B33')) + geom_jitter()

jpeg('SCFA_levels_barplots.jpg',units = "in", width = 5.2, height = 4, res = 300)
gridExtra::grid.arrange(acet, prop, butyr, nrow = 1, ncol = 3)
dev.off()
