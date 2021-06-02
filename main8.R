# Libraries

# Functions
getLOH <- function(tab) {
  df <- as.data.frame(table(tab))
  num <- df[df$Var1 == "D",]$Freq + df[df$Var1 == "L",]$Freq
  loh <- 100*num/sum(df$Freq)
  return(loh)
}

# Get the file to draw the plots from the command line
sysarg <- commandArgs(trailingOnly = TRUE)
file <- sysarg[1][1]

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

#, length(pos$submitter), " positive, ", length(neg$submitter), " negative, and ", length(neu$submitter), " neutral\n")
