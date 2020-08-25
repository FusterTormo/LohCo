# Constantes. Nombres esperados de los archivos para hacer los graficos y los calculos de coverage
lista <- '../gensAestudi.txt'
base <- 'coverageBase.txt'
all <- 'coverageAll.txt'
salida <- 'coverage.txt' # Archivo donde se guardan las estadisticas de coverage

#Crear un grafico con el porcentaje de bases que tienen un coverage X. En el eje X se dibuja el coverage y en el eje Y el porcentaje de bases que tienen dicho coverage
coverageGeneral <- function(cov){
  #Calcular coverage maximo a dibujar a partir del coverage medio
  maxCov <- round(max(cov$coverage)*1.5)

  #Leer los datos
  cov <- read.table(all)
  colnames(cov) <- c("etiqueta","coverage","bases amb aquest coverage","total bases","% de bases amb aquest coverage")
  #Sumar las fracciones de coverage que tiene cada region
  gcov_cumul = 1 - cumsum(cov[,5])
  png(filename="coverage.png",width=800,height=800)
  #Dibujar el coverage
  plot(cov[2:(maxCov+1),2], gcov_cumul[1:maxCov], type='l',xlab='Coverage',ylab=paste("%reads","amb","coverage"),col='dodgerblue',ylim=c(0,1),lwd=1,axes=FALSE)
  #Dibujar los ejes de coordenadas
  aux <- par("usr")
  axis(2, at=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'), las=1)
  rang <- seq(from=0, to=aux[2], by=200)
  axis(1, at=rang, labels=rang)
  box()
  #Lineas adicionales para que el grafico sea mas legible
  for(it in rang) {
    abline(v=it,col="gray60")
  }
  for(it in seq(from=0,to=1,by=0.1)){
    abline(h=it,col="gray60")
  }
  dev.off()
}

# Dibuja un grafico con el coverage que tiene cada uno de los genes que hay en el archivo gensAestudi.txt. Asume que existen, en el directorio de trabajo, los archivos gensAestudi.txt y coveragePerBase.txt
coveragePerGen <- function(cov){
  # Leer la lista de genes. Se crea un grafico por cada fila
  gens <- read.table(lista)
  colnames(gens) <- c("nom")
  # Lista con los coverages minimo, maximo, medio y mediana
  stats <- data.frame()
  # Crear los graficos. Uno por cada gen que hay en la lista de genes (gensAestudi.txt)
  for(it in 1:length(gens$nom)) {
    aux <- cov[grep(gens$nom[it],cov$gen),]
    maxCov <- max(aux$coverage)
    minCov <- min(aux$coverage)
    meanCov <- mean(aux$coverage)
    medianCov <- median(aux$coverage)
    aux2 <- data.frame(gene=gens$nom[it], min=minCov, max=maxCov, mean=meanCov, median=medianCov)
    stats <- rbind(stats,aux2)
    png(filename=paste("coverage",gens$nom[it],".png",sep=""),width=500,height=600)
    plot(aux$coverage, type="l", xlab="Bases", ylab="Coverage", ylim=c(0,maxCov), col="deepskyblue1", lwd=2, main=paste("Coverage", "de", gens$nom[it]))
    abline(h=0.0, col = "gray20")
    abline(h=100, col = "firebrick1") # Rojo
    abline(h=500, col = "darkorange") # Naranja
    abline(h=1000, col = "forestgreen") # Verde
    dev.off()
  }
  # Escribir la lista con los coverages minimo, maximo, medio y mediana en un archivo de texto
  write.table(stats, sep = "\t", file = "coverageGeneStats.tsv", row.names = FALSE)
}

# Extraer el minimo, maximo, la media, la mediana de coverage, así como el porcentaje de bases cubiertos a 0, 30, 100, 500 y 1000 de coverage. Los datos se guardan en un archivo de texto llamado coverage.txt
estadisticas <- function(cov) {
  minimo <- min(cov$coverage)
  maximo <- max(cov$coverage)
  media <- mean(cov$coverage)
  mediana <- median(cov$coverage)
  total <- length(cov$coverage)
  cov0 <- length(cov[cov$coverage == 0,]$coverage)
  cov30 <- length(cov[cov$coverage <= 30,]$coverage)
  cov100 <- length(cov[cov$coverage <= 100,]$coverage)
  cov500 <- length(cov[cov$coverage <= 500,]$coverage)
  cov1000 <- length(cov[cov$coverage <= 1000,]$coverage)
  cat("{\'minimo\':", minimo, ",", file = salida)
  cat("\'maximo\':", maximo, ",", file = salida, append = TRUE)
  cat("\'media\':", media, ",", file = salida, append = TRUE)
  cat("\'mediana\':", mediana, ",", file = salida, append = TRUE)
  cat("\'bases0\':", cov0/total*100, ",", file = salida, append = TRUE)
  cat("\'bases30\':", cov30/total*100, ",", file = salida, append = TRUE)
  cat("\'bases100\':", cov100/total*100, ",", file = salida, append = TRUE)
  cat("\'bases500\':", cov500/total*100, ",", file = salida, append = TRUE)
  cat("\'bases1000\':", cov1000/total*100, "}\n", file = salida, append = TRUE)
}

#Programa principal
cov <- read.table(base)
colnames(cov) <- c("chr","start","end","gen","base","coverage")
estadisticas(cov)
coverageGeneral(cov)
coveragePerGen(cov)
