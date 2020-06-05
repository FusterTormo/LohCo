library(karyoploteR)
wd.temp <- getwd()
setwd("../Desktop/comparacionsFVS/images")
wd <- "../"

createPlots <- function() {
  files <- list.files(wd, full.names = TRUE, recursive = FALSE, pattern = "*txt")
  lapply(files, function(path) {
    name <- strsplit(basename(path), ".txt")[1]
    txt <- read.csv(path, header = FALSE)
    same <- as.numeric(txt[1,])
    diff <- as.numeric(txt[2,])
    #aux <- paste("raw", name, sep = "_")
    #figname <- paste(aux, "png", sep = ".")
    #png(figname, width = 900, height = 900)
    #boxplot(same, diff, names = c("Same", "Different"), xlab = "Regions", ylab = "Region size", main = name)
    #dev.off()
    aux <- paste("box", name, sep = "_")
    figname <- paste(aux, "png", sep = ".")
    png(figname, width = 900, height = 900)
    boxplot(same, diff, outline = FALSE, names = c("Same", "Different"), col = c("green", "red"), xlab = "Regions", ylab = "Region size", main = name)
    dev.off()
    aux <- paste("line", name, sep = "_")
    figname <- paste(aux, "png", sep = ".")
    png(figname, width = 900, height = 900)
    plot(sort(same), pch = 18, col = "green", xlab = "Regions", ylab = "Region size", main = name)
    points(sort(diff), col = "red", pch = 18)
    dev.off()
  })
}

checkMedians <- function() {
  # Comprobar si la mediana de la longitud de las regiones que coinciden es menor que la mediana de longitud de las regiones que no coinciden
  files <- list.files(wd, full.names = TRUE, recursive = FALSE, pattern = "*txt")
  fs <- read.table("../Desktop/facetsVSsequenza.tsv", sep = "\t", header = TRUE)
  lapply(files, function(path) {
    name <- basename(path)
    txt <- read.csv(path, header = FALSE)
    same <- as.numeric(txt[1,])
    diff <- as.numeric(txt[2,])
    s <- as.numeric(summary(same)[3]) # Get the median for the regions that have the same aberration
    d <- as.numeric(summary(diff)[3])
    if (s < d) {
      case <- strsplit(name, ".txt")[1]
      print(fs[fs$Case == case, c("Case", "regSim", "baseSim")])
    }
  })
}

plotKaryotypes <- function() {
  
  files <- list.files(wd, full.names = TRUE, recursive = FALSE, pattern = "*regsCoin.bed")
  lapply(files, function(coin) {
    case <- strsplit(basename(coin), "regsCoin.bed")[1]
    diff <- gsub("regsCoin", "regsDiff", coin)
    coin.df <- read.table(coin, header = FALSE, sep = "\t")
    colnames(coin.df) <- c("chr", "start", "end")
    coin.df$chr <- paste0("chr", coin.df$chr)
    coin.gr <- toGRanges(coin.df)
    diff.df <- read.table(diff, header = FALSE, sep = "\t")
    colnames(diff.df) <- c("chr", "start", "end")
    diff.df$chr <- paste0("chr", diff.df$chr)
    diff.gr <- toGRanges(diff.df)
    aux <- paste("karyo", case, sep = "_")
    figname <- paste(aux, "png", sep = ".")
    png(figname, width = 900, height = 900)
    kp <- plotKaryotype(chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
    kpPlotRegions(kp, data = coin.gr, col = "#00FF00", layer.margin = 0.01, border = NA, r0 = 0.1, r1 = 0.3)
    kpPlotRegions(kp, data = diff.gr, col = "#FF0000", layer.margin = 0.01, border = NA, r0 = 0.4, r1 = 0.6)
    dev.off()
  })
}

plotKaryotypes()
createPlots()
setwd(wd.temp)