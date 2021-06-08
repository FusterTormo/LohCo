#################################################################
# Checking the combination purity and mean copy number reported #
#################################################################
library(ggplot2)
library(scales)

# Load the data
full <- read.table("meanCN.tsv", sep = "\t", header = TRUE)

# Purity correlation between LOH tools
purity <- data.frame(FACETS = full$fac_purity, Sequenza = full$seq_purity, PURPLE = full$pur_purity, ascatNGS = full$ngs_purity)
ggplot(purity, aes(FACETS, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggplot(purity, aes(FACETS, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggplot(purity, aes(FACETS, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggplot(purity, aes(Sequenza, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggplot(purity, aes(Sequenza, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggplot(purity, aes(PURPLE, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()

# Mean CN correlation
cn <- data.frame(ASCAT2 = full$asc_meanCN, FACETS = full$fac_meanCN, Sequenza = full$seq_meanCN, PURPLE = full$pur_meanCN, ascatNGS = full$ngs_meanCN)
ggplot(cn, aes(ASCAT2, FACETS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(ASCAT2, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(ASCAT2, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(ASCAT2, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(FACETS, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(FACETS, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(FACETS, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(Sequenza, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(Sequenza, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggplot(cn, aes(PURPLE, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()

# Copy number reported
boxplot(cn, col = hue_pal()(5), main = "Mean copy number reported")
boxplot(cn, col = hue_pal()(5), main = "Mean copy number reported", outline = FALSE)

# Aberration percentage reported by each tool
colors <- hue_pal()(4)
tmp <- strsplit(full$asc_aberration, ";")
ascat <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
boxplot(ascat, names = c("Amp", "LOH", "Del", "WT"), col = colors, main = "ASCAT2")
tmp <- strsplit(full$fac_aberration, ";")
facets <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
boxplot(facets, names = c("Amp", "LOH", "Del", "WT"), col = colors, main = "FACETS")
tmp <- strsplit(full$seq_aberration, ";")
sequenza <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
boxplot(sequenza, names = c("Amp", "LOH", "Del", "WT"), col = colors, main = "Sequenza")
tmp <- strsplit(full$pur_aberration, ";")
purple <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
boxplot(purple, names = c("Amp", "LOH", "Del", "WT"), col = colors, main = "PURPLE")
tmp <- strsplit(full$ngs_aberration, ";")
ascatngs <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
boxplot(ascatngs, names = c("Amp", "LOH", "Del", "WT"), col = colors, main = "ascatNGS")

