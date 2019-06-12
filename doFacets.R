showHelp <- function() {
  cat("Usage:\n\n")
  cat("\tRscript doFacets.R folder_with_FACETS folder [reference_genome] [cval_for_preProcSample cval_for_procSample]\n")
  cat("WHERE\n\t- reference_genome is the reference genome, which can be \"hg19\", \"hg38\", \"hg18\", \"mm9\", \"mm10\", \"udef\"\n")
  cat("\t- cval_preProcSample is the value of cval for preProcSample function\n")
  cat("\t- cval_procSample is the value of cval for procSample function\n\n")
}

executeFacets <- function(folder, geno, cval1, cval2) {
  library ("facets")
  filename <- "facets_comp.csv"
  basicFile <- "facets_comp_basic.tsv"
  cncfFile <- "facets_comp_cncf.tsv"
  pngFile <- "facets_comp.png"
  spiderPlot <- "facets_spider.png"
  cat("INFO: Executing FACETS in ", filename , " \n")
  set.seed(1234)
  rcmat = readSnpMatrix(filename)
  xx = preProcSample(rcmat, gbuild = geno, cval = cval1)
  oo = procSample(xx,cval = cval2)
  fit = emcncf(oo)
  cat("INFO: Storing info")
  write.table(fit$cncf, file=cncfFile, sep="\t", row.names = FALSE)
  basic <- data.frame(fit$loglik, fit$purity, fit$ploidy, fit$dipLogR)
  colnames(basic) <- c("Log_likelihood", "Purity", "Ploidy", "Estimated_logR_diploid_segments")
  write.table(basic, file=basicFile, sep="\t", row.names = FALSE)
  png(pngFile, width = 900, height = 900, units = "px", pointsize = 18)
  plotSample(oo, emfit = fit)
  dev.off()
  png(spiderPlot)
  logRlogORspider(oo$out, oo$dipLogR)
  dev.off()
  cat ("INFO: Data stored succesfully:\n\tfacets_comp_cncf -> Contains the columns of segmentation output\n\tfacets_comp_basic.tsv -> Contains information about purity, ploidy ...\n")
  cat ("\tfacets_comp.png -> Contains the plot generated using FACETS plotSample command\n")
}

#Get the arguments from the command line
sysarg <- commandArgs(trailingOnly = TRUE)
#Constants
gen <- "hg38"
cval1 <- 25
cval2 <- 150

if (length(sysarg) == 0 || length(sysarg) > 4) {
  showHelp()
}else if(length(sysarg) == 4) {
  dir <- sysarg[1]
  gen <- sysarg[2]
  cval1 <- as.numeric(sysarg[3])
  cval2 <- as.numeric(sysarg[4])
}else if (length(sysarg) == 2) {
  dir <- sysarg[1]
  cval1 <- as.numeric(sysarg[2])
  cval2 <- as.numeric(sysarg[3])
}else if(length(sysarg) == 2) {
  dir <- sysarg[1]
  gen <- sysarg[2]
}else {
  dir <- sysarg[1]
}

if (!is.na(cval1) && !is.na(cval2)) {
  executeFacets(dir, gen, cval1, cval2)
}