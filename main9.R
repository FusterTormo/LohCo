library(knitr)

# Functions
printSummary <- function(data) {
  cat("R-INFO: Summary stats:\n")
  cat("\tMean: ", mean(data), " Median: ", median(data), " Min: ", min(data), " Max: ", max(data), "\n")
}

# Main program
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
gr <- read.table("grouped.vars.tsv", header = TRUE, sep = "\t")
cat("R-INFO: Data loaded successfully. Creating the plots\n")

# Assert the coordinates in the histograms are the same
# First position
if (p$position[1] != t$position[1])
  cat("R-WARNING: Initial coordinates in positive and pathogenic are not the same\n")
if (p$position[1] != n$position[1])
  cat("R-WARNING: Initial coordinates in positive and negative are not the same\n")
if (t$position[1] != n$position[1])
  cat("R-WARNING: Initial coordinates in negative and pathogenic are not the same\n")
# Last position
if (rev(p$position)[1] != rev(t$position)[1])
  cat("R-WARNING: Final coordinates in positive and pathogenic are not the same\n")
if (rev(p$position)[1] != rev(n$position)[1])
  cat("R-WARNING: Final coordinates in positive and negative are not the same\n")
if (rev(t$position)[1] != rev(n$position)[1])
  cat("R-WARNING: Final coordinates in negative and pathogenic are not the same\n")

# Plots
cat("R-INFO: Creating the plots\n")
# Times a position is mutated in the gene
png(width = 1261, height = 906, filename = paste(gene, "_Pos.png", sep = ""))
plot(p$position, p$times, pch = 18, col = "red", main = "Variant position", xlab = paste(gene, "position"), ylab = "Variants found", xlim = c(min(p$position, n$position, t$position), max(p$position, n$position, t$position)))
points(n$position, n$times, pch = 18, col = "blue")
points(t$position, t$times, pch = 18, col = "orange")
legend("topleft", fill = c("red", "blue", "orange"), legend = c("Positive", "Negative", "Pathogenic"))
dev.off()
# Frequency that a position is mutated in the gene
png(width = 1261, height = 906, filename = paste(gene, "_Freq.png", sep = ""))
plot(p$position, p$freq, pch = 18, col = "red", main = "Variant position", xlab = paste(gene, "position"), ylab = "Variant frequency", xlim = c(min(p$position, n$position, t$position), max(p$position, n$position, t$position)))
points(n$position, n$freq, pch = 18, col = "blue")
points(t$position, t$freq, pch = 18, col = "orange")
legend("topleft", fill = c("red", "blue", "orange"), legend = c("Positive", "Negative", "Pathogenic"))
dev.off()
# Select the positions in positive submitters where there are variants
pzoom <- p[p$freq > 0,]
tzoom <- t[t$position %in% pzoom$position,]
nzoom <- n[n$position %in% pzoom$position,]
png(width = 1261, height = 906, filename = paste(gene, "_ZoomFreq.png", sep = ""))
plot(pzoom$freq, pch = 18, col = "red", main = "Variant position", xlab = "BRCA1 position", ylab = "Variants found", ylim = c(0, max(pzoom$freq)))
points(tzoom$freq, pch = 18, col = "orange")
points(nzoom$freq, pch = 18, col = "blue")
legend("topleft", fill = c("red", "blue", "orange"), legend = c("Positive", "Negative", "Pathogenic"))
dev.off()

# Barplot with variant type
tp <- as.table(table(vp$type))
png(width=1261, height = 906, filename = paste(gene, "_VarsPositive.png", sep = ""))
barplot(100*prop.table(tp), col = rainbow(length(tp)), main = "Variant type percent in positive submitters", ylim = c(0, 100))
dev.off()
tt <- as.table(table(vt$type))
png(width=1261, height = 906, filename = paste(gene, "_VarsPathogenic.png", sep = ""))
barplot(100*prop.table(tt), col = rainbow(length(tt)), main = "Variant type percent in pathogenic submitters", ylim = c(0, 100))
dev.off()
tn <- as.table(table(vn$type))
png(width=1261, height = 906, filename = paste(gene, "_VarsNegative.png", sep = ""))
barplot(100*prop.table(tn), col = rainbow(length(tn)), main = "Variant type percent in negative submitters", ylim = c(0, 100))
dev.off()

# Variants per submitter
png(width=1261, height = 941, filename = paste(gene, "_varsXsubmitter.png", sep = ""))
boxplot(as.matrix(table(vp$submitter))[,1], as.matrix(table(vt$submitter))[,1], as.matrix(table(vn$submitter))[,1], col = c("red", "orange", "blue"), main = "Variants per submitter", names = c("Positive", "Pathogenic", "Negative"))
dev.off()

# Clinvar reported significance
png(paste(gene, "_significance.png"), width = 1261, height = 941)
barplot(table(vp$cln.signf), col = "red", main = "Reported significance in positive group")
barplot(table(vt$cln.signf), col = "orange", main = "Reported significance in pathogenic group")
barplot(table(vn$cln.signf), col = "blue", main = "Reported significance in negative group")
dev.off()

# Tables
# Variants present in LOH positive/pathogenic submitters but not in LOH negative submitters
d <- gr[gr$InLOHNegative == 0,]
if (length(d$Chr) < 5) {
  cat("R-INFO: Variants present in LOH positive/pathogenic submitters, but not in LOH negative submitters\n")
  knitr::kable(d) # Print the table in console, in a fancy manner
}
write.table(d, file = "allCandVars.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

e <- d[d$Type == "exonic" | d$Type == "splicing",]
if (length(e$Chr) < 5) {
  cat("\n\nR-INFO: Coding variants present in LOH positive/pathogenic submitters, but not in LOH negative submitters\n")
  knitr::kable(e) # Print the table in console, in a fancy manner
}
write.table(e, file = "candVars.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cat("\nR-INFO: Variants in ", gene, " reported for each LOH positive submitter\n")
printSummary(as.matrix(table(vp$submitter))[,1])

cat("R-INFO: Variants in ", gene, " reported for each LOH pathogenic submitter\n")
printSummary(as.matrix(table(vt$submitter))[,1])

cat("R-INFO: Variants in ", gene, " reported for each LOH negative submitter\n")
printSummary(as.matrix(table(vn$submitter))[,1])

cat("R-INFO: Script finished successfully. Check plots, specially ZoomFreq, for further clues. \n\tCandidate variants stored as allCandVars.tsv and candVars.tsv. \n\tRun main9.py in this folder to search this variants in other cancer repositories\n")

#####  Old code ##### 
# # Get the number of variants per submitter
# e1 <- as.data.frame(table(pos$submitter))
# e2 <- as.data.frame(table(neg$submitter))
# colnames(e1) <- c("submitter", "variants")
# colnames(e2) <- c("submitter", "variants")
# ggp <- data.frame(name = "Positive", variants = e1$variants)
# ggn <- data.frame(name = "Negative", variants = e2$variants)
# gg <- rbind(ggp, ggn)
# ggplot(gg, aes(y=variants, x=name, fill = name)) + geom_violin(trim = FALSE) + geom_boxplot(width = 0.1) + theme_classic() + scale_x_discrete(name = "") + ggtitle("Variants per submitter")
# ggsave(filename = paste(gene, "_varsXsubmitter.png", sep = ""), device = "png", dpi = 320, width = 26, height = 26, units = "cm")
# 
# 

# # Descriptive statistics of the number of variants per submitter
# tmp <- as.data.frame(table(pos$submitter))
# cat("R-INFO: Mean variants per submitter in positive cases: ", mean(tmp$Freq))
# cat("R-INFO: Median variants per submitter in positive cases: ", median(tmp$Freq))
# tmp <- as.data.frame(table(neg$submitter))
# cat("R-INFO: Mean variants per submitter in positive cases: ", mean(tmp$Freq))
# cat("R-INFO: Median variants per submitter in positive cases: ", median(tmp$Freq))
# 
# # Calculate dN/dS ratio in positive/negative samples
# non <- exonic[exonic$name == "nonsynonymous SNV",][2:3]
# syn <- exonic[exonic$name == "synonymous SNV",][2:3]
# aux <- non/syn
# cat("R-INFO: dN/dS ratio for positive cases: ", non$positive , "/", syn$positive, "= ", aux$positive, "\n")
# cat("R-INFO: dN/dS ratio for negative cases: ", non$negative , "/", syn$negative, "= ", aux$negative, "\n")
