#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
MAIN: Anotar el fichero de coverages creado por GATK. Usar dicho fichero para crear las graficas de coverage de la forma en que se hace usando bedtools + coveragePanells.R
"""

import os
import sys

import constantes as cte

def guardaManifest(mani) :
    """Guardar el manifest en un diccionario y devolverlo

    Parameters
    ----------
        mani : str
            Ruta del manifest que se quiere anotar. Se supone que no lleva cabecera

    Returns
    -------
        dict
            Diccionario con el formato: {"crom" : [{"inicio" : int, "fin" : int, "gen" : str}]}
    """
    ret = {}
    with open(mani, "r") as fi :
        for l in fi :
            aux = l.strip().split("\t")
            if len(aux) > 0 :
                crom = aux[0]
                ini = int(aux[1])
                fin = int(aux[2])
                gen = aux[3]
                if crom in ret.keys() :
                    ret[crom].append({"inicio" : ini, "fin" : fin, "gen" : gen})
                else :
                    ret[crom] = [{"inicio" : ini, "fin" : fin, "gen" : gen}]
    return ret

def getGene(crom, pos, mani) :
    """Buscar en el manifest el gen que corresponde la posicion y el cromosoma pasado por parametro

    Parameters
    ----------
        crom : str
            Nombre del cromosoma. Debe tener el mismo formato que las claves del diccionario del manifest
        pos : int
            Posicion genomica a buscar en el manifest
        mani : dict
            Manifest guardado con el formato generado usando guardarManifest

    Returns
    -------
        str
            Nombre del gen que se ha encontrado. En caso de no encontrar la posicion dentro del manifest, devolvera ""
    """
    gen = ""
    if crom in mani.keys() :
        pos = int(pos)
        for l in mani[crom] :
            if pos >= l["inicio"] and pos <= l["fin"] :
                gen = l["gen"]
                break
    else :
        print("WARNING: {} no encontrado en {}".format(crom, mani.keys()))

    return gen

def anotar(arx, mani = cte.manifest) :
    """Anotar el archivo de GATK

    Anota los genes de cada posicion en el archivo de coverage creado por GATK DepthOfCoverage. Usa un manifest pasado por parametro, o el manifest que hay en la libreria
    de constantes. Crea 2 archivos:
        1. Con el mismo contenido que tenia GATK mas una columna adicional con el gen al que pertenece la posicion genomica (columna Locus)
        2. Un archivo con los datos agrupados por coverage. Es decir, cuantas bases tienen coverage=500, cuantas 501, etc.

    Parameters
    ----------
        arx : str
            Ruta donde esta el archivo GATK DepthOfCoverage que se quiere anotar. Las columnas deben estar separadas por tabulaciones
        mani : str, optional
            Ruta del manifest a partir del cual se anotaran las coordenadas

    Returns
    -------
        str
            Nombre del archivo con los coverages anotados
        str
            Nombre del archivo con las bases agrupadas por coverage
        list
            Lista de genes en los que se ha buscado coverage
    """
    manifest = {} # Contenido del manifest
    cnt = "" # Contenido que tendra el archivo de coverage anotado
    genes = []
    coverages = {}
    bases = 0
    output = ""
    if os.path.isfile(mani) :
        print("INFO: Leyendo manifest")
        manifest = guardaManifest(mani)
        header = "" # Se supone que el archivo de coverage tiene cabecera
        if os.path.isfile(arx) :
            print("INFO: Leyendo archivo de coverage")
            output = "{}.anno.cov".format(arx)
            output2 = "{}.all.cov".format(arx)
            with open(arx, "r") as fi :
                for l in fi :
                    if header == "" :
                        header = "Locus\tDepth\tAvg_depth\tDepth_sample\tBase_counts\tGene\n"
                    else :
                        aux = l.strip().split("\t")
                        crom, pos = aux[0].split(":")
                        tmp = int(aux[1])
                        bases += 1
                        if tmp not in coverages.keys() :
                            coverages[tmp] = 1
                        else :
                            coverages[tmp] += 1
                        gen = getGene(crom, pos, manifest)
                        if gen not in genes :
                            genes.append(gen)
                        cnt += "{}\t{}\n".format(l.strip(), gen)
            print("INFO: Guardando los datos como {}".format(output))
            with open(output, "w") as fi :
                fi.write(header)
                fi.write(cnt)
            print("INFO: Y como {}".format(output2))
            keys = coverages.keys()
            with open(output2, "w") as fi :
                fi.write("tag\tcoverage\tbases\ttotal_bases\tpercent\n")
                for k in sorted(keys) :
                    fi.write("all\t{cov}\t{times}\t{all}\t{div}\n".format(cov = k, times = coverages[k], all = bases, div = coverages[k]/bases))
        else :
            print("ERROR: Archivo de coverage no encontrado en {}".format(arx))
            sys.exit(1)
    else :
        print("ERROR: Manifest no encontrado en {}".format(mani))
        sys.exit(1)

    return output, output2, genes

def crearScript(arx, arx2, genes) :
    """Crear un archivo de R que haga los graficos de coverage a partir de los datos de los archivos pasados por parametro

    Parameters
    ----------
        arx : str
            Ruta del archivo con formato GATK DepthOfCoverage anotado. Este parametro lo genera la funcion "anotar"
        arx2 : str
            Ruta del archivo con las bases agrupadas por coverage. Este parametro lo genera la funcion "anotar"
        genes : list
            Lista de genes en los que dibujar el coverage individual
    """
    print("INFO: Creado script R")
    script = "plotCoverage.R"
    with open(script, "w") as fi :
        fi.write("#Leer los datos de coverage\n")
        fi.write("all <- read.table('{}', sep = '\\t', header = TRUE)\n".format(arx))
        fi.write("general <- read.table('{}', sep = '\\t', header = TRUE)\n".format(arx2))
        fi.write("\n#-------Calcular estadisticas descriptivas-------\n")
        fi.write("minimo <- min(all$Depth)\n")
        fi.write("maximo <- max(all$Depth)\n")
        fi.write("media <- mean(all$Depth)\n")
        fi.write("mediana <- median(all$Depth)\n")
        fi.write("total <- length(all$Depth)\n")
        fi.write("cov0 <- length(all[all$Depth == 0,]$Depth)\n")
        fi.write("cov30 <- length(all[all$Depth <= 30,]$Depth)\n")
        fi.write("cov100 <- length(all[all$Depth <= 100,]$Depth)\n")
        fi.write("cov500 <- length(all[all$Depth <= 500,]$Depth)\n")
        fi.write("cov1000 <- length(all[all$Depth <= 1000,]$Depth)\n")
        fi.write("cat(\"{\'minimo\':\", minimo, \",\", file = '../coverage.txt')\n")
        fi.write("cat(\"\'maximo\':\", maximo, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'media\':\", media, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'mediana\':\", mediana, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'bases0\':\", cov0/total*100, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'bases30\':\", cov30/total*100, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'bases100\':\", cov100/total*100, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'bases500\':\", cov500/total*100, \",\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("cat(\"\'bases1000\':\", cov1000/total*100, \"}\\n\", file = '../coverage.txt', append = TRUE)\n")
        fi.write("\n#-------Coverage general-------\n")
        fi.write("maxCov <- max(general$coverage)\n")
        fi.write("gcov.cumul <- 1 - cumsum(general$percent)\n")
        fi.write("png('coverage.png', width = 800, height = 800)\n")
        fi.write("plot(general[2:(maxCov+1),2], gcov.cumul[1:maxCov], type = 'l', xlab='Coverage',ylab=paste('%reads','amb','coverage'),col='dodgerblue', ylim=c(0,1), lwd=1, axes=FALSE, main = 'Coverage')\n")
        fi.write("# Dibujar los ejes de coordenadas\n")
        fi.write("aux <- par('usr')\n")
        fi.write("axis(2, at=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'), las=1)\n")
        fi.write("rang <- seq(from=0, to=aux[2], by=200)\n")
        fi.write("axis(1, at=rang, labels=rang)\n")
        fi.write("box()\n")
        fi.write("# Lineas adicionales para que sea mas legible\n")
        fi.write("for(it in rang) {\n")
        fi.write("\tabline(v=it,col='gray60')\n")
        fi.write("}\n")
        fi.write("for(it in seq(from=0,to=1,by=0.1)) {\n")
        fi.write("\tabline(h=it,col='gray60')\n")
        fi.write("}\n")
        fi.write("dev.off()\n")
        for g in genes :
            fi.write("\n#-------Coverage {}-------\n".format(g))
            fi.write("temp <- all[all$Gene == '{}',]\n".format(g))
            fi.write("maxCov <- max(temp$Depth)\n")
            fi.write("png('coverage{}.png', width = 500, height = 600)\n".format(g))
            fi.write("plot(temp$Depth, type = 'l', xlab = 'Bases', ylab = 'Coverage', ylim = c(0, maxCov), col = 'deepskyblue1', lwd=2, main=paste('Coverage', 'de', '{}'))\n".format(g))
            fi.write("abline(h=0.0, col = 'gray20')\n")
            fi.write("abline(h=100, col = 'firebrick1') # Rojo\n")
            fi.write("abline(h=500, col = 'darkorange') # Naranja\n")
            fi.write("abline(h=1000, col = 'forestgreen') # Verde\n")
            fi.write("dev.off()\n")

"""Programa principal"""
if __name__ == "__main__" :
    if len(sys.argv) > 1 :
        if len(sys.argv) > 2 :
            archivo, archivo2, genes = anotar(sys.argv[1], sys.argv[2])
        else :
            archivo, archivo2, genes = anotar(sys.argv[1])
        crearScript(archivo, archivo2, genes)
    else :
        arx = input("INPUT: Sobre que archivo anotar el coverage? ")
        aux = input("INPUT: Usar este manifest {} ? (S/N) ".format(cte.manifest))
        if aux.lower() == "s" :
            archivo, archivo2, genes = anotar(arx)
        else :
            mani = input("INPUT: Ruta del manifest a usar para anotar: ")
            archivo, archivo2, genes = anotar(arx, mani)
        crearScript(archivo, archivo2, genes)
