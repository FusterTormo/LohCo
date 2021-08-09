library(ggplot2)
# Get the gene name for 
sysarg <- commandArgs(trailingOnly = TRUE)
gene = sysarg[1][1]

# Load the data
cat("R-INFO: Loading data\n")
p <- read.table("positiveHistogram.tsv", header = TRUE, sep = "\t")
t <- read.table("pathogenicHistogram.tsv", header = TRUE, sep = "\t")
n <- read.table("negativeHistogram.tsv", header = TRUE, sep = "\t")
vp <- read.table("posVariants.tsv", header = TRUE, sep = "\t")
vt <- read.table("patVariants.tsv", header = TRUE, sep = "\t")
vn <- read.table("negVariants.tsv", header = TRUE, sep = "\t")
cat("R-INFO: Data loaded successfully. Creating the plots\n")

# Assert the coordinates in the histograms are the same
# First position
if (p$position[1] != t$position[1])
  cat("R-WARNING: Initial coordinates in positive and pathogenic are not the same")
if (p$position[1] != n$position[1])
  cat("R-WARNING: Innitial coordinates in positive and negative are not the same")
if (t$position[1] != n$position[1])
  cat("R-WARNING: Initial coordinates in negative and pathogenic are not the same")
# Last position
if (rev(p$position)[1] != rev(t$position)[1])
  cat("R-WARNING: Final coordinates in positive and pathogenic are not the same")
if (rev(p$position)[1] != rev(n$position)[1])
  cat("R-WARNING: Final coordinates in positive and negative are not the same")
if (rev(t$position)[1] != rev(n$position)[1])
  cat("R-WARNING: Final coordinates in negative and pathogenic are not the same")

# Plots

# Number of times a position is mutated in the gene
png(width = 1261, height = 906, filename = paste(gene, "_genePos.png", sep = ""))
plot(p$times, pch = 18, col = "red", main = "Variant position", xlab = paste(gene, "position"), ylab = "Variants found")
points(n$times, pch = 18, col = "blue")
points(t$times, pch = 18, col = "orange")
legend("topleft", fill = c("red", "blue", "orange"), legend = c("Positive", "Negative", "Pathogenic"))
dev.off()
# Frequency that this position is mutated in the gene
png(width = 1261, height = 906, filename = paste(gene, "_geneFreq.png", sep = ""))
plot(p$freq, pch = 18, col = "red", main = "Variant position", xlab = paste(gene, "position"), ylab = "Variant frequency")
points(n$freq, pch = 18, col = "blue")
points(t$freq, pch = 18, col = "orange")
legend("topleft", fill = c("red", "blue", "orange"), legend = c("Positive", "Negative", "Pathogenic"))
dev.off()

# Barplot with variant type
tab <- table(vp$type)
total <- length(vp$type)
perc <- 100*tab/total
barplot(perc, col = rainbow(length(tab)), main = "Variant type percent in positive submitters")
tab <- table(vt$type)
total <- length(vt$type)
perc <- 100*tab/total
barplot(perc, col = rainbow(length(tab)), main = "Variant type percent in pathogenic submitters")
tab <- table(vn$type)
total <- length(vn$type)
perc <- 100*tab/total
barplot(perc, col = rainbow(length(tab)), main = "Variant type percent in negative submitters")

# Variants per submitter
boxplot(table(vp$submitter), table(vt$submitter), table(vn$submitter), col = c("red", "orange", "blue"), main = "Variants per submitter")

# Clinvar reported significance
barplot(table(vp$cln.signf), col = "red", main = "Reported significance in positive group")
barplot(table(vt$cln.signf), col = "orange", main = "Reported significance in pathogenic group")
barplot(table(vn$cln.signf), col = "blue", main = "Reported significance in negative group")

#TODO continuar per aci per si queda alguna grafica per transferir
# TODO editar els colors de les grafiques i fer les comandes per guardar en el disc
# The plots


# Zoom by selecting only the positions with more than 1 variant
pzoom <- p[p$times > 1,]
nzoom <- n[n$position %in% pzoom$position,]
png(width = 1261, height = 906, filename = paste(gene, "_zoom_genePos.png", sep = ""))
plot(pzoom$times, pch = 18, col = "red", main = "Variant position", xlab = "BRCA1 position", ylab = "Variants found")
points(nzoom$times, pch = 20, col = "blue")
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

# Plot the SNV changes
# TODO. Create the matrix with all the C>A, A>G... changes
# TODO. Do a beside barplot with this matrix
snv <- paste0(pos$ref, ">", pos$alt)
posSNV <- as.data.frame(table(snv))
snv <- paste0(neg$ref, ">", neg$alt)
negSNV <- as.data.frame(table(snv))

# Type of variants
# TODO. Get the exonic variant type per group (positive, negative)

# Descriptive statistics of the number of variants per submitter
tmp <- as.data.frame(table(pos$submitter))
cat("R-INFO: Mean variants per submitter in positive cases: ", mean(tmp$Freq))
cat("R-INFO: Median variants per submitter in positive cases: ", median(tmp$Freq))
tmp <- as.data.frame(table(neg$submitter))
cat("R-INFO: Mean variants per submitter in positive cases: ", mean(tmp$Freq))
cat("R-INFO: Median variants per submitter in positive cases: ", median(tmp$Freq))

# Calculate dN/dS ratio in positive/negative samples
non <- exonic[exonic$name == "nonsynonymous SNV",][2:3]
syn <- exonic[exonic$name == "synonymous SNV",][2:3]
aux <- non/syn
cat("R-INFO: dN/dS ratio for positive cases: ", non$positive , "/", syn$positive, "= ", aux$positive, "\n")
cat("R-INFO: dN/dS ratio for negative cases: ", non$negative , "/", syn$negative, "= ", aux$negative, "\n")
