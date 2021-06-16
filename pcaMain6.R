library(ggbiplot)

# Do PCA regarding JCC obtained using ASCAT2 as truth set
data <- read.table("ascat2_heatmap.tsv", sep = "\t", header = TRUE)
data.clean <- data[!is.na(data$ascat2) & !is.na(data$dnacopy) & !is.na(data$facets) & !is.na(data$ascatngs) & !is.na(data$sequenza) & !is.na(data$purple),]
mt <- matrix(c(data.clean$ascat2, data.clean$dnacopy, data.clean$facets, data.clean$ascatngs, data.clean$sequenza, data.clean$purple), ncol = 6)
rownames(mt) <- data.clean$submitter
colnames(mt) <- c("ASCAT2", "DNAcopy", "FACETS", "ascatNGS", "Sequenza", "PURPLE")
# Plot heatmap
png("ascat2_heatmap.png", width = 922, height = 620)
heatmap(mt, main = "Jaccard index using ASCAT2 as truth set")
dev.off()
# Do tool PCA
# First remove ASCAT2 column as it has no variance and raises problems when standardizing the data
mt <- matrix(c(data.clean$dnacopy, data.clean$facets, data.clean$ascatngs, data.clean$sequenza, data.clean$purple), ncol = 5)
rownames(mt) <- data.clean$submitter
colnames(mt) <- c("DNAcopy", "FACETS", "ascatNGS", "Sequenza", "PURPLE")
png("ascat2_heatmap2.png", width = 922, height = 620)
heatmap(mt, main = "Jaccard index using ASCAT2 as truth set")
dev.off()
mt.pca <- prcomp(mt, center = TRUE, scale. = TRUE)
ggbiplot(mt.pca, obs.scale = TRUE) + theme_minimal() + ggtitle("Jaccard Index compared with ASCAT2")
ggsave("ascat2_pca.png")
#ggbiplot(mt.pca, obs.scale = T, labels = rownames(mt))

# Do PCA regarding JCC obtained using ASCAT2 as truth set
data <- read.table("dnacopy_heatmap.tsv", sep = "\t", header = TRUE)
data.clean <- data[!is.na(data$ascat2) & !is.na(data$dnacopy) & !is.na(data$facets) & !is.na(data$ascatngs) & !is.na(data$sequenza) & !is.na(data$purple),]
mt <- matrix(c(data.clean$ascat2, data.clean$dnacopy, data.clean$facets, data.clean$ascatngs, data.clean$sequenza, data.clean$purple), ncol = 6)
rownames(mt) <- data.clean$submitter
colnames(mt) <- c("ASCAT2", "DNAcopy", "FACETS", "ascatNGS", "Sequenza", "PURPLE")
# Plot heatmap
png("dnacopy_heatmap.png", width = 922, height = 620)
heatmap(mt, main = "Jaccard index using DNAcopy as truth set")
dev.off()
# Do tool PCA
# First remove DNAcopy column as it has no variance and raises problems when standardizing the data
mt <- matrix(c(data.clean$ascat2, data.clean$facets, data.clean$ascatngs, data.clean$sequenza, data.clean$purple), ncol = 5)
rownames(mt) <- data.clean$submitter
colnames(mt) <- c("ASCAT2", "FACETS", "ascatNGS", "Sequenza", "PURPLE")
png("dnacopy_heatmap2.png", width = 922, height = 620) 
heatmap(mt, main = "Jaccard index using DNAcopy as truth set")
dev.off()
mt.pca <- prcomp(mt, center = TRUE, scale. = TRUE)
ggbiplot(mt.pca, obs.scale = TRUE) + theme_minimal() + ggtitle("Jaccard Index compared with DNAcopy")
ggsave("dnacopy_pca.png")