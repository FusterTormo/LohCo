# Libraries

# Functions
getLOH <- function(tab) {
  df <- as.data.frame(table(tab, useNA = "always"))
  dels <- df[df$tab == "D" & !is.na(df$tab),]$Freq
  loss <- df[df$tab == "L" & !is.na(df$tab),]$Freq
  if (length(dels) <= 0)
    dels <- 0
  if (length(loss) <= 0)
    loss <- 0
  num <- dels  + loss
  loh <- round(100*num/sum(df$Freq),2)
  return(loh)
}

convert2vector <- function(tab) {
  df <- as.data.frame(table(tab))
  a <- df[df$tab == "A",]$Freq
  if (length(a) <= 0)
    a <- 0
  l <- df[df$tab == "L",]$Freq
  if (length(l) <= 0)
    l <- 0
  d <- df[df$tab == "D",]$Freq
  if (length(d) <= 0)
    d <- 0
  n <- df[df$tab == "N",]$Freq
  if (length(n) <= 0)
    n <- 0
  nf <- length(tab) - a - l - d - n
  vec <- c(a, d, l, n, nf)
  return(vec)
}

# Get the file to draw the plots from the command line
sysarg <- commandArgs(trailingOnly = TRUE)
file <- sysarg[1][1]

# Get the gene name. It would be printed in the input file as the first characters before "_"
gene <- strsplit(file, "_")[[1]][1]

# Load data
cat("R-INFO: Loading data\n")
all <- read.table(file, header = TRUE, sep = "\t")
pos <- all[all$GermlineVar == '+',]
neg <- all[all$GermlineVar == "-",]
neu <- all[all$GermlineVar == "?",]

cat("R-INFO: ", length(all$submitter), " total cases\n")

cat("R-INFO: ", length(pos$submitter), " considered positive\n")
cat("\tASCAT2 reported ", getLOH(pos$ASCAT2), "% LOH\n")
cat("\tFACETS reported ", getLOH(pos$FACETS), "% LOH\n")
cat("\tascatNGS reported ", getLOH(pos$ascatNGS), "% LOH\n")
cat("\tPURPLE reported ", getLOH(pos$PURPLE), "% LOH\n")
cat("\tSequenza reported ", getLOH(pos$Sequenza), "% LOH\n")
cat("\tLOH custom reported ", getLOH(pos$LOHcat), "% LOH\n")

cat("R-INFO: ", length(neg$submitter), " considered negative\n")
cat("\tASCAT2 reported ", getLOH(neg$ASCAT2), "% LOH\n")
cat("\tFACETS reported ", getLOH(neg$FACETS), "% LOH\n")
cat("\tascatNGS reported ", getLOH(neg$ascatNGS), "% LOH\n")
cat("\tPURPLE reported ", getLOH(neg$PURPLE), "% LOH\n")
cat("\tSequenza reported ", getLOH(neg$Sequenza), "% LOH\n")
cat("\tLOH custom reported ", getLOH(neg$LOHcat), "% LOH\n")

cat("R-INFO: ", length(pos$submitter), " considered neutral\n")
cat("\tASCAT2 reported ", getLOH(neu$ASCAT2), "% LOH\n")
cat("\tFACETS reported ", getLOH(neu$FACETS), "% LOH\n")
cat("\tascatNGS reported ", getLOH(neu$ascatNGS), "% LOH\n")
cat("\tPURPLE reported ", getLOH(neu$PURPLE), "% LOH\n")
cat("\tSequenza reported ", getLOH(neu$Sequenza), "% LOH\n")
cat("\tLOH custom reported ", getLOH(neu$LOHcat), "% LOH\n")

cat("R-INFO: Without the variant classification\n")
cat("\tASCAT2 reported ", getLOH(all$ASCAT2), "% LOH\n")
cat("\tFACETS reported ", getLOH(all$FACETS), "% LOH\n")
cat("\tascatNGS reported ", getLOH(all$ascatNGS), "% LOH\n")
cat("\tPURPLE reported ", getLOH(all$PURPLE), "% LOH\n")
cat("\tSequenza reported ", getLOH(all$Sequenza), "% LOH\n")
cat("\tLOH custom reported ", getLOH(all$LOHcat), "% LOH\n")

# Convert the data to vectors
pa <- convert2vector(pos$ASCAT2)
pf <- convert2vector(pos$FACETS)
pn <- convert2vector(pos$ascatNGS)
pp <- convert2vector(pos$PURPLE)
ps <- convert2vector(pos$Sequenza)
pl <- convert2vector(pos$LOHcat)

na <- convert2vector(neg$ASCAT2)
nf <- convert2vector(neg$FACETS)
nn <- convert2vector(neg$ascatNGS)
np <- convert2vector(neg$PURPLE)
ns <- convert2vector(neg$Sequenza)
nl <- convert2vector(neg$LOHcat)

ua <- convert2vector(neu$ASCAT2)
uf <- convert2vector(neu$FACETS)
un <- convert2vector(neu$ascatNGS)
up <- convert2vector(neu$PURPLE)
us <- convert2vector(neu$Sequenza)
ul <- convert2vector(neu$LOHcat)

# Convert the data to matrix (one matrix per tool to create the barplots)
ma <- matrix(c(pa, na, ua), nrow = 3, byrow = TRUE)
mf <- matrix(c(pf, nf, uf), nrow = 3, byrow = TRUE)
mn <- matrix(c(pn, nn, un), nrow = 3, byrow = TRUE)
mp <- matrix(c(pp, np, up), nrow = 3, byrow = TRUE)
ms <- matrix(c(ps, ns, us), nrow = 3, byrow = TRUE)
ml <- matrix(c(pl, nl, ul), nrow = 3, byrow = TRUE)

# Constants that will be used in the barplots
top <- max(ma, mf, mn, mp, ms, ml) # Maximum for the ylim parameter
colors <- c("red", "green", "yellow") # Color for barplots
xtag <- c("Amp", "LOH", "Del", "Norm", "Not_Av") # Tags for the names

# Create barplots with the absolute data
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = paste("ASCAT2 aberrations reported in", gene), ylim = c(0,top))
barplot(mf, beside = TRUE, col = colors, names.arg = xtag, main = paste("FACETS aberrations reported in", gene), ylim = c(0,top))
barplot(mn, beside = TRUE, col = colors, names.arg = xtag, main = paste("ascatNGS aberrations reported in", gene), ylim = c(0,top))
barplot(mp, beside = TRUE, col = colors, names.arg = xtag, main = paste("PURPLE aberrations reported in", gene), ylim = c(0,top))
barplot(ms, beside = TRUE, col = colors, names.arg = xtag, main = paste("Sequenza aberrations reported in", gene), ylim = c(0,top))
barplot(ml, beside = TRUE, col = colors, names.arg = xtag, main = paste("Custom LOH aberrations reported in", gene), ylim = c(0,top))

# Convert the data to percentages (normalize the data)
ma[1,] <- ma[1,]/sum(ma[1,])*100
ma[2,] <- ma[2,]/sum(ma[2,])*100
ma[3,] <- ma[3,]/sum(ma[3,])*100
mf[1,] <- mf[1,]/sum(mf[1,])*100
mf[2,] <- mf[2,]/sum(mf[2,])*100
mf[3,] <- mf[3,]/sum(mf[3,])*100
mn[1,] <- mn[1,]/sum(mn[1,])*100
mn[2,] <- mn[2,]/sum(mn[2,])*100
mn[3,] <- mn[3,]/sum(mn[3,])*100
mp[1,] <- mp[1,]/sum(mp[1,])*100
mp[2,] <- mp[2,]/sum(mp[2,])*100
mp[3,] <- mp[3,]/sum(mp[3,])*100
ms[1,] <- ms[1,]/sum(ms[1,])*100
ms[2,] <- ms[2,]/sum(ms[2,])*100
ms[3,] <- ms[3,]/sum(ms[3,])*100
ml[1,] <- ml[1,]/sum(ml[1,])*100
ml[2,] <- ml[2,]/sum(ml[2,])*100
ml[3,] <- ml[3,]/sum(ml[3,])*100

# Create the barplots with the normalized data
barplot(ma, beside = TRUE, col = colors, names.arg = xtag, main = paste("ASCAT2 aberrations reported in", gene), ylim = c(0,100))
barplot(mf, beside = TRUE, col = colors, names.arg = xtag, main = paste("FACETS aberrations reported in", gene), ylim = c(0,100))
barplot(mn, beside = TRUE, col = colors, names.arg = xtag, main = paste("ascatNGS aberrations reported in", gene), ylim = c(0,100))
barplot(mp, beside = TRUE, col = colors, names.arg = xtag, main = paste("PURPLE aberrations reported in", gene), ylim = c(0,100))
barplot(ms, beside = TRUE, col = colors, names.arg = xtag, main = paste("Sequenza aberrations reported in", gene), ylim = c(0,100))
barplot(ml, beside = TRUE, col = colors, names.arg = xtag, main = paste("Custom LOH aberrations reported in", gene), ylim = c(0,100))
