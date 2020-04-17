coverageGeneral <- function(resum){
  #Calcular coverage maximo a dibujar a partir del coverage medio
  maxCov <- round(resum[4]*1.5)
  #Leer los datos
  cov <- read.table('coverageAll.txt')
  colnames(cov) <- c("etiqueta","coverage","bases amb aquest coverage","total bases","% de bases amb aquest coverage")
  #Sumar las fracciones de coverage que tiene cada region
  gcov_cumul = 1 - cumsum(cov[,5])
  png(filename="coverage.png",width=800,height=800)
  #Dibujar el coverage (maximo 1000)
  plot(cov[2:(maxCov+1),2], gcov_cumul[1:maxCov], type='l',xlab='Coverage',ylab=paste("%reads","amb","coverage"),col='blue',ylim=c(0,1),lwd=3,axes=FALSE)
  #Dibujar los ejes de coordenadas
  aux <- par("usr")
  axis(2,at=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'),las=1)
  rang <- seq(from=0,to=aux[2],by=200)
  axis(1,at=rang,labels=rang)
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

coveragePerGen <- function(){
  gens <- read.table('gensAestudi.txt')
  colnames(gens) <- c("nom")
  cov <- read.table('coveragePerBase.txt')
  if (length(cov) == 6){
    colnames(cov) <- c("chr","start","end","gen","base","coverage")
  } else {
    colnames(cov) <- c("chr","start","end","gen","alt","strand","base","coverage")
  }
  for(it in 1:length(gens$nom)) {
    png(filename=paste("coverage",gens$nom[it],".png",sep=""),width=500,height=600)
    aux <- cov[grep(gens$nom[it],cov$gen),]
    maxCov <- summary(aux$coverage)[6]
    plot(aux$coverage, type="l",xlab="Bases",ylab="Coverage",ylim=c(0,maxCov),col="blue",main=paste("Coverage","de",gens$nom[it]),lwd=2)
    abline(h=0.0,col="gray60")
    abline(h=200,col="gray60")
    dev.off()
  }
}

percentatges <- function(){
  cov <- read.table('coveragePerBase.txt')
  if (length(cov) == 6){
    colnames(cov) <- c("chr","start","end","gen","base","coverage")
  } else {
    colnames(cov) <- c("chr","start","end","gen","alt","strand","base","coverage")
  }
  covTotal <- length(cov$end)
  covZero <- length(cov[cov$coverage == 0,]$end)
  cov30 <- length(cov[cov$coverage <= 30,]$end)
  cov100 <- length(cov[cov$coverage <= 100,]$end)
  cov500 <- length(cov[cov$coverage <= 500,]$end)
  covHigh <- length(cov[cov$coverage <= 1000,]$end)
  dades <- c("Bases_coverage_=_0",covZero,"Percentatge",(covZero/covTotal)*100,"Bases_coverage_<=_30",cov30,"Percentatge",(cov30/covTotal)*100,"Bases_coverage_<=_100",cov100,"Percentatge",(cov100/covTotal)*100,"Bases_coverage_<=_500",cov500,"Percentatge",(cov500/covTotal)*100,"Bases_coverage_<=_1000",covHigh,"Percentatge",(covHigh/covTotal)*100,"Bases_totals",covTotal)
  resum <- summary(cov$coverage)
  cat("Coverage mig:", resum[4], "X\n")
  cat("Coverage minim:", resum[1],"\n")
  cat("Coverage maxim:", resum[6],"\n")
  write(dades,file="estadistiquesCoverage.txt",ncolumns=4)
  return (resum)
}

#Programa principal
coverages <- percentatges()
coverageGeneral(coverages)
coveragePerGen()
