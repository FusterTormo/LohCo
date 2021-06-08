##############################################################################################
# MAIN: Plot the comparisons between ASCAT2, ascatNGS, DNAcopy, FACETS, PURPLE, and Sequenza #
##############################################################################################

library(ggplot2)
library(vioplot)
library(scales)

# Codes
# d -> DNAcopy
# a -> ASCAT2
# n -> ascatNGS
# f -> FACETS
# s -> Sequenza
# p -> PURPLE

colors <- hue_pal()(6)
setwd("main6/")

# Load data
# DNAcopy vs all
dva <- read.table("arrayVSascat2.tsv", header = TRUE, sep = "\t")
dvd <- read.table("arrayVSarray.tsv", header = TRUE, sep = "\t")
dvf <- read.table("arrayVSfacets.tsv", header = TRUE, sep = "\t")
dvn <- read.table("arrayVSascatNGS.tsv", header = TRUE, sep = "\t")
dvs <- read.table("arrayVSsequenza.tsv", header = TRUE, sep = "\t")
dvp <- read.table("arrayVSpurple.tsv", header = TRUE, sep = "\t")
# ASCAT2 vs all
ava <- read.table("ascat2VSascat2.tsv", header = TRUE, sep = "\t")
avd <- read.table("ascat2VSarray.tsv", header = TRUE, sep = "\t")
avf <- read.table("ascat2VSfacets.tsv", header = TRUE, sep = "\t")
avn <- read.table("ascat2VSascatNGS.tsv", header = TRUE, sep = "\t")
avs <- read.table("ascat2VSsequenza.tsv", header = TRUE, sep = "\t")
avp <- read.table("ascat2VSpurple.tsv", header = TRUE, sep = "\t")
# FACETS vs all
fva <- read.table("facetsVSascat2.tsv", header = TRUE, sep = "\t")
fvd <- read.table("facetsVSarrays.tsv", header = TRUE, sep = "\t")
fvf <- read.table("facetsVSfacets.tsv", header = TRUE, sep = "\t")
fvn <- read.table("facetsVSascatNGS.tsv", header = TRUE, sep = "\t")
fvs <- read.table("facetsVSsequenza.tsv", header = TRUE, sep = "\t")
fvp <- read.table("facetsVSpurple.tsv", header = TRUE, sep = "\t")
# ascatNGS vs all
nva <- read.table("ascatNGSVSascat2.tsv", header = TRUE, sep = "\t")
nvd <- read.table("ascatNGSVSarrays.tsv", header = TRUE, sep = "\t")
nvf <- read.table("ascatNGSVSfacets.tsv", header = TRUE, sep = "\t")
nvn <- read.table("ascatNGSVSascatNGS.tsv", header = TRUE, sep = "\t")
nvs <- read.table("ascatNGSVSsequenza.tsv", header = TRUE, sep = "\t")
nvp <- read.table("ascatNGSVSpurple.tsv", header = TRUE, sep = "\t")
# Sequenza vs all
sva <- read.table("sequenzaVSascat2.tsv", header = TRUE, sep = "\t")
svd <- read.table("sequenzaVSarrays.tsv", header = TRUE, sep = "\t")
svf <- read.table("sequenzaVSfacets.tsv", header = TRUE, sep = "\t")
svn <- read.table("sequenzaVSascatNGS.tsv", header = TRUE, sep = "\t")
svs <- read.table("sequenzaVSsequenza.tsv", header = TRUE, sep = "\t")
svp <- read.table("sequenzaVSpurple.tsv", header = TRUE, sep = "\t")
# PURPLE vs all
pva <- read.table("purpleVSascat2.tsv", header = TRUE, sep = "\t")
pvd <- read.table("purpleVSarrays.tsv", header = TRUE, sep = "\t")
pvf <- read.table("purpleVSfacets.tsv", header = TRUE, sep = "\t")
pvn <- read.table("purpleVSascatNGS.tsv", header = TRUE, sep = "\t")
pvs <- read.table("purpleVSsequenza.tsv", header = TRUE, sep = "\t")
pvp <- read.table("purpleVSpurple.tsv", header = TRUE, sep = "\t")

# Violin plots of whole Jaccard Index
datos <- rbind(data.frame(JCC = dvd$JCC, grp = "DNAcopy"), data.frame(JCC = dva$JCC, grp = "ASCAT2"), 
               data.frame(JCC = dvf$JCC, grp = "FACETS"), data.frame(JCC = dvn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = dvs$JCC, grp = "Sequenza"), data.frame(JCC = dvp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("DNAcopy vs other tools")
ggsave("DNAcopyVSall.png")

datos <- rbind(data.frame(JCC = avd$JCC, grp = "DNAcopy"), data.frame(JCC = ava$JCC, grp = "ASCAT2"), 
               data.frame(JCC = avf$JCC, grp = "FACETS"), data.frame(JCC = avn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = avs$JCC, grp = "Sequenza"), data.frame(JCC = avp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("ASCAT2 vs other tools")
ggsave("ASCAT2VSall.png")

datos <- rbind(data.frame(JCC = fvd$JCC, grp = "DNAcopy"), data.frame(JCC = fva$JCC, grp = "ASCAT2"), 
               data.frame(JCC = fvf$JCC, grp = "FACETS"), data.frame(JCC = fvn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = fvs$JCC, grp = "Sequenza"), data.frame(JCC = fvp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("FACETS vs other tools")
ggsave("FACETSVSall.png")

datos <- rbind(data.frame(JCC = nvd$JCC, grp = "DNAcopy"), data.frame(JCC = nva$JCC, grp = "ASCAT2"), 
               data.frame(JCC = nvf$JCC, grp = "FACETS"), data.frame(JCC = nvn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = nvs$JCC, grp = "Sequenza"), data.frame(JCC = nvp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("ascatNGS vs other tools")
ggsave("ascatNGSVSall.png")

datos <- rbind(data.frame(JCC = svd$JCC, grp = "DNAcopy"), data.frame(JCC = sva$JCC, grp = "ASCAT2"), 
               data.frame(JCC = svf$JCC, grp = "FACETS"), data.frame(JCC = svn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = svs$JCC, grp = "Sequenza"), data.frame(JCC = svp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("Sequenza vs other tools")
ggsave("SequenzaVSall.png")

datos <- rbind(data.frame(JCC = pvd$JCC, grp = "DNAcopy"), data.frame(JCC = pva$JCC, grp = "ASCAT2"), 
               data.frame(JCC = pvf$JCC, grp = "FACETS"), data.frame(JCC = pvn$JCC, grp = "ascatNGS"), 
               data.frame(JCC = pvs$JCC, grp = "Sequenza"), data.frame(JCC = pvp$JCC, grp = "PURPLE"))
ggplot(datos, aes(grp, JCC, fill = grp)) + geom_violin() + geom_boxplot(width = 0.1) + theme_minimal() + 
  theme(legend.position = "none") + expand_limits(y = c(-0.1, 1.1)) + xlab("") + ylab("Jaccard Index") + 
  ggtitle("PURPLE vs other tools")
ggsave("PURPLEVSall.png")

# Violin plots of Jaccard Index by aberration
# Amplifications
#plot(density(dva$MCCamp), type = "l", col = colors[1])
#lines(density(dvn$MCCamp), col = colors[2])
# Violin plots of Matthew's Correlation Coefficient by aberration