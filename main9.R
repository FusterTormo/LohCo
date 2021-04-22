library(ggplot2)
# Get the gene name for 
sysarg <- commandArgs(trailingOnly = TRUE)
gene = sysarg[1][1]

# Load the data
cat("R-INFO: Loading data\n")
p <- read.table("positionHistogram.tsv", header = TRUE, sep = "\t")
n <- read.table("negativeHistogram.tsv", header = TRUE, sep = "\t")
mt <- matrix(c(p$times, n$times), nrow = 2, byrow = TRUE)
pos <- read.table("posVariants.tsv", header = FALSE, sep = "\t")
neg <- read.table("negVariants.tsv", header = FALSE, sep = "\t")
colnames(pos) <- c("chr", "start", "end", "ref", "alt", "position", "exonic", "submitter", "tm", "cn")
colnames(neg) <- c("chr", "start", "end", "ref", "alt", "position", "exonic", "submitter", "tm", "cn")

cat("R-INFO: Data loaded successfully. Creating the plots\n")
# The plots
# Number of times a position is mutated in BRCA1
png(width = 1261, height = 906, filename = paste(gene, "_genePos.png", sep = ""))
plot(p$times, pch = 18, col = "red", main = "Variant position", xlab = "BRCA1 position", ylab = "Variants found")
points(n$times, pch = 20, col = "blue")
legend("topleft", fill = c("red", "blue"), legend = c("Positive", "Negative"))
dev.off()

# Barplot different mutations (exonic, intronic, UTR...)
png(width = 1261, height = 906, filename = paste(gene, "_varType.png", sep = ""))
e1 <- as.data.frame(table(pos$position))
e2 <- as.data.frame(table(neg$position))
both <- merge(e1, e2, by = "Var1", all = TRUE)
both[is.na(both)] <- 0 # Change 'NA' to 0
mt <- t(matrix(c(both$Freq.x, both$Freq.y), ncol = 2))
p <- barplot(mt, beside = TRUE, main = "Variants position", col = c("red", "blue"), names.arg = both$Var1, legend.text = c("Positive", "Negative"), ylim = c(0, max(mt)+100))
text(p, mt+40, labels = mt)
dev.off()

# Barplot different exonic mutations
png(width = 1261, height = 906, filename = paste(gene, "_exonicVarType.png", sep = ""))
e1 <- as.data.frame(table(pos$exonic))
e2 <- as.data.frame(table(neg$exonic))
exonic <- merge(e1, e2, by = "Var1", all = TRUE)
exonic[is.na(exonic)] <- 0  # Change 'NA' to 0
colnames(exonic) <- c("name", "positive", "negative")
mt <- t(matrix(c(exonic$positive, exonic$negative), ncol = 2))
p <- barplot(mt, beside = TRUE, main = "Exonic variants", col = c("red", "blue"), names.arg = exonic$name, legend.text = c("Positive", "Negative"), ylim = c(0, max(mt)+30))
text(p, mt+10, labels = mt)
dev.off()

# Get the number of variants per submitter
e1 <- as.data.frame(table(pos$submitter))
e2 <- as.data.frame(table(neg$submitter))
colnames(e1) <- c("submitter", "variants")
colnames(e2) <- c("submitter", "variants")
ggp <- data.frame(name = "Positive", variants = e1$variants)
ggn <- data.frame(name = "Negative", variants = e2$variants)
gg <- rbind(ggp, ggn)
ggplot(gg, aes(y=variants, x=name, fill = name)) + geom_violin(trim = FALSE) + geom_boxplot(width = 0.1) + theme_classic() + scale_x_discrete(name = "") + ggtitle("Variants per submitter")
ggsave(filename = paste(gene, "_varsXsubmitter.png", sep = ""), device = "png", dpi = 320, width = 26, height = 26, units = "cm")

# Calculate dN/dS ratio in positive/negative samples
non <- exonic[exonic$name == "nonsynonymous SNV",][2:3]
syn <- exonic[exonic$name == "synonymous SNV",][2:3]
aux <- non/syn
cat("R-INFO: dN/dS ratio for positive cases: ", non$positive , "/", syn$positive, "= ", aux$positive, "\n")
cat("R-INFO: dN/dS ratio for positive cases: ", non$negative , "/", syn$negative, "= ", aux$negative, "\n")