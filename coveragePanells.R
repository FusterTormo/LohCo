#Crear un grafico con el porcentaje de bases que tienen un coverage X. En el eje X se dibuja el coverage y en el eje Y el porcentaje de bases que tienen dicho coverage
coverageGeneral <- function(){
  # Esta funcion asume que, en el directorio de trabajo, existen los archivos coveragePerBase.txt y coverageAll.txt
  base <- 'coveragePerBase.txt'
  all <- 'coverageAll.txt'
  
  #Calcular coverage maximo a dibujar a partir del coverage medio
  cov <- read.table(base)
  colnames(cov) <- c("chr","start","end","gen","base","coverage")
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
coveragePerGen <- function(){
  lista <- 'gensAestudi.txt'
  base <- 'coveragePerBase.txt'
  # Leer la lista de genes. Se crea un grafico por cada fila
  gens <- read.table(lista)
  colnames(gens) <- c("nom")
  # Leer el archivo de coverage por base
  cov <- read.table(base)
  colnames(cov) <- c("chr","start","end","gen","base","coverage")
  
  # Crear los graficos. Uno por cada gen que hay en la lista de genes (gensAestudi.txt)
  for(it in 1:length(gens$nom)) {
    png(filename=paste("coverage",gens$nom[it],".png",sep=""),width=500,height=600)
    aux <- cov[grep(gens$nom[it],cov$gen),]
    maxCov <- max(aux$coverage)
    plot(aux$coverage, type="l", xlab="Bases", ylab="Coverage", ylim=c(0,maxCov), col="deepskyblue1", lwd=2, main=paste("Coverage", "de", gens$nom[it]))
    abline(h=0.0, col = "gray20")
    abline(h=100, col = "gray30")
    abline(h=500, col = "gray40")
    abline(h=1000, col = "gray50")
    dev.off()
  }
}

#Programa principal
coverageGeneral()
coveragePerGen()
