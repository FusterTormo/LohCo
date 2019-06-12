accuracy <- function (ruta) {
  a <- read.table(ruta, sep = "\t", header = TRUE)
  acc <- (a$Del_FA[1]+a$Norm_FA[2]+a$Amp_FA[3])/(sum(a$Del_FA)+sum(a$Norm_FA)+sum(a$Amp_FA))
  return (acc)
}

file.list <- list.files("./", full.names = TRUE, recursive = FALSE)
vc <- c()
for (i in 1:length(file.list)) {
  aux <- accuracy(file.list[i])
  vc <- c(vc, aux)
}

cat (summary(vc))

