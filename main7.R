#################################################################
# Checking the combination purity and mean copy number reported #
#################################################################
library(ggplot2)
library(scales)

# Load the data
full <- read.table("meanCN.tsv", sep = "\t", header = TRUE)

# Purity correlation between LOH tools
auxPurity <- data.frame(FACETS = full$fac_purity, Sequenza = full$seq_purity, PURPLE = full$pur_purity, ascatNGS = full$ngs_purity)
purity <- auxPurity[!is.na(auxPurity$FACETS) & !is.na(auxPurity$Sequenza) & !is.na(auxPurity$PURPLE) & !is.na(auxPurity$ascatNGS),]
ggplot(purity, aes(FACETS, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("purityFvS.png")
ggplot(purity, aes(FACETS, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("purityFvP.png")
ggplot(purity, aes(FACETS, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("purityFvN.png")
ggplot(purity, aes(Sequenza, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("puritySvP.png")
ggplot(purity, aes(Sequenza, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("puritySvN.png")
ggplot(purity, aes(PURPLE, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Purity correlation") + theme_minimal()
ggsave("purityPvN.png")

# Mean CN vs ploidy correlation
auxPl <- data.frame(ploidy = full$fac_ploidy, meancn = full$fac_meanCN)
fac <- auxPl[!is.na(auxPl$ploidy) & !is.na(auxPl$meancn),]
ggplot(fac, aes(ploidy, meancn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() + ggtitle("FACETS ploidy vs mean CN") + xlab("Ploidy") + ylab("Mean CN")
ggsave("FACETSploVScn.png")
auxPl <- data.frame(ploidy = full$seq_ploidy, meancn = full$seq_meanCN)
seq <- auxPl[!is.na(auxPl$ploidy) & !is.na(auxPl$meancn),]
ggplot(seq, aes(ploidy, meancn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() + ggtitle("Sequenza ploidy vs mean CN") + xlab("Ploidy") + ylab("Mean CN")
ggsave("SEQUENZAploVScn.png")
auxPl <- data.frame(ploidy = full$ngs_ploidy, meancn = full$ngs_meanCN)
ngs <- auxPl[!is.na(auxPl$ploidy) & !is.na(auxPl$meancn),]
ggplot(ngs, aes(ploidy, meancn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() + ggtitle("ascatNGS ploidy vs mean CN") + xlab("Ploidy") + ylab("Mean CN")
ggsave("ASCATNGSploVScn.png")
auxPl <- data.frame(ploidy = full$pur_ploidy, meancn = full$pur_meanCN)
pur <- auxPl[!is.na(auxPl$ploidy) & !is.na(auxPl$meancn),]
ggplot(pur, aes(ploidy, meancn)) + geom_point() + geom_smooth(method = "lm") + theme_minimal() + ggtitle("PURPLE ploidy vs mean CN") + xlab("Ploidy") + ylab("Mean CN")
ggsave("PURPLEploVScn.png")

# Mean CN correlation
auxCn <- data.frame(ASCAT2 = full$asc_meanCN, FACETS = full$fac_meanCN, Sequenza = full$seq_meanCN, PURPLE = full$pur_meanCN, ascatNGS = full$ngs_meanCN)
cn <- auxCn[!is.na(auxCn$ASCAT2) & !is.na(auxCn$FACETS) & !is.na(auxCn$Sequenza) & !is.na(auxCn$PURPLE) & !is.na(auxCn$ascatNGS),]
ggplot(cn, aes(ASCAT2, FACETS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnAvF.png")
ggplot(cn, aes(ASCAT2, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnAvS.png")
ggplot(cn, aes(ASCAT2, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnAvP.png")
ggplot(cn, aes(ASCAT2, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnAvN.png")
ggplot(cn, aes(FACETS, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnFvS.png")
ggplot(cn, aes(FACETS, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnFvP.png")
ggplot(cn, aes(FACETS, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnFvN.png")
ggplot(cn, aes(Sequenza, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnSvP.png")
ggplot(cn, aes(Sequenza, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnSvN.png")
ggplot(cn, aes(PURPLE, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Mean copy number") + theme_minimal()
ggsave("cnPvS.png")

# Ploidy correlation
auxPloidy <- data.frame(FACETS = full$fac_ploidy, Sequenza = full$seq_ploidy, PURPLE = full$pur_ploidy, ascatNGS = full$ngs_ploidy)
ploidy <- auxPloidy[!is.na(auxPloidy$FACETS) & !is.na(auxPloidy$Sequenza) & !is.na(auxPloidy$PURPLE) & !is.na(auxPloidy$ascatNGS),]
ggplot(ploidy, aes(FACETS, Sequenza)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plFvS.png")
ggplot(ploidy, aes(FACETS, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plFvP.png")
ggplot(ploidy, aes(FACETS, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plFvN.png")
ggplot(ploidy, aes(Sequenza, PURPLE)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plSvP.png")
ggplot(ploidy, aes(Sequenza, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plSvN.png")
ggplot(ploidy, aes(PURPLE, ascatNGS)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Ploidy") + theme_minimal()
ggsave("plPvS.png")

# Copy number reported. Notice that data for boxplot comes from the no NAs data
png("allCN.png", width = 720, height = 554)
boxplot(cn, col = hue_pal()(5), main = "Mean copy number reported")
dev.off()
png("allCN_outliers.png", width = 720, height = 554)
boxplot(cn, col = hue_pal()(5), main = "Mean copy number reported", outline = FALSE)
dev.off()

# Ploidy reported. Same happens here than before
png("ploidyCN.png", width = 720, height = 554)
boxplot(ploidy, col = hue_pal()(4), main = "Ploidy reported")
dev.off()

# Aberration percentage reported by each tool
colors <- hue_pal()(4)
tmp.nas <- full[full$asc_aberration != "NA;NA;NA;NA",]$asc_aberration
tmp <- strsplit(as.character(tmp.nas), ";")
mt <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
ascat <- rbind(data.frame(AB = "Amplification", percent = mt[,1]), 
               data.frame(AB = "LOH", percent = mt[,2]),
               data.frame(AB = "Deletion", percent = mt[,3]),
               data.frame(AB = "Wild-type", percent = mt[,4]))
ggplot(ascat, aes(AB, percent, fill = AB)) + geom_violin() + theme_minimal() + theme(legend.position = "none") + ggtitle("ASCAT2") + xlab("Aberration") + ylab("Aberration percent per sample")
ggsave("ascat2Aberrations.png")

tmp.nas <- full[full$fac_aberration != "NA;NA;NA;NA",]$fac_aberration
tmp <- strsplit(as.character(tmp.nas), ";")
mt <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
facets <- rbind(data.frame(AB = "Amplification", percent = mt[,1]), 
               data.frame(AB = "LOH", percent = mt[,2]),
               data.frame(AB = "Deletion", percent = mt[,3]),
               data.frame(AB = "Wild-type", percent = mt[,4]))
ggplot(facets, aes(AB, percent, fill = AB)) + geom_violin() + theme_minimal() + theme(legend.position = "none") + ggtitle("FACETS") + xlab("Aberration") + ylab("Aberration percet per sample")
ggsave("facetsAberrations.png")

tmp.nas <- full[full$seq_aberration != "NA;NA;NA;NA",]$seq_aberration
tmp <- strsplit(as.character(tmp.nas), ";")
mt <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
sequenza <- rbind(data.frame(AB = "Amplification", percent = mt[,1]), 
               data.frame(AB = "LOH", percent = mt[,2]),
               data.frame(AB = "Deletion", percent = mt[,3]),
               data.frame(AB = "Wild-type", percent = mt[,4]))
ggplot(sequenza, aes(AB, percent, fill = AB)) + geom_violin() + theme_minimal() + theme(legend.position = "none") + ggtitle("Sequenza") + xlab("Aberration") + ylab("Aberration percet per sample")
ggsave("sequenzaAberrations.png")

tmp.nas <- full[full$pur_aberration != "NA;NA;NA;NA",]$pur_aberration
tmp <- strsplit(as.character(tmp.nas), ";")
mt <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
purple <- rbind(data.frame(AB = "Amplification", percent = mt[,1]), 
               data.frame(AB = "LOH", percent = mt[,2]),
               data.frame(AB = "Deletion", percent = mt[,3]),
               data.frame(AB = "Wild-type", percent = mt[,4]))
ggplot(purple, aes(AB, percent, fill = AB)) + geom_violin() + theme_minimal() + theme(legend.position = "none") + ggtitle("PURPLE") + xlab("Aberration") + ylab("Aberration percet per sample")
ggsave("purpleAberrations.png")

tmp.nas <- full[full$ngs_aberration != "NA;NA;NA;NA",]$ngs_aberration
tmp <- strsplit(as.character(tmp.nas), ";")
mt <- matrix(as.numeric(unlist(tmp)), ncol = 4, byrow = TRUE)
ascatngs <- rbind(data.frame(AB = "Amplification", percent = mt[,1]), 
               data.frame(AB = "LOH", percent = mt[,2]),
               data.frame(AB = "Deletion", percent = mt[,3]),
               data.frame(AB = "Wild-type", percent = mt[,4]))
ggplot(ascatngs, aes(AB, percent, fill = AB)) + geom_violin() + theme_minimal() + theme(legend.position = "none") + ggtitle("ascatNGS") + xlab("Aberration") + ylab("Aberration percet per sample")
ggsave("ascatngsAberrations.png")

