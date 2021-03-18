#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys

import constantes as cte
import vcfQC

def recogerDatos(ruta = "./") :
    if ruta == "./" :
        print("INFO: Realizando control de calidad en {}".format(os.getcwd()))
    else :
        print("INFO: Realizando control de calidad en {}".format(ruta))
    # Leer las carpetas con los nombres de las muestras
    for root, dirs, files in os.walk(ruta) :
        break

    print("INFO: {} muestras encontradas".format(len(dirs)))
    datos = {}
    for d in dirs :
        datos[d] = {"FQ" : 0, "BAM" : 0, "ON" : 0, "OFF" : 0, "COV" : {}, "VCF" : {}, "COV_PATH" : ""}
        path = "{}/{}".format(ruta, d)
        qc = "{}/{}".format(path, cte.qcaln)
        # Guardar en el diccionario los parametros de calidad del bam
        if os.path.isfile(qc) :
            with open(qc, "r") as fi :
                aux = fi.read()
            aux2 = eval(aux)
            datos[d]["FQ"] = aux2["FASTQ"]
            datos[d]["BAM"] = aux2["BAM"]
            datos[d]["ON"] = aux2["ON"]
            datos[d]["OFF"] = aux2["OFF"]
        else :
            print("WARNING: {} no encontrado. Ejecuta AUP/bamQC.py para generarlo".format(qc))

        # Guardar en el diccionario los parametros de calidad del coverage
        cov = "{}/{}".format(path, cte.covarx)
        if os.path.isfile(cov) :
            with open(cov, "r") as fi :
                aux = fi.read()
            aux2 = eval(aux)
            datos[d]["COV"]["min"] = aux2["minimo"]
            datos[d]["COV"]["max"] = aux2["maximo"]
            datos[d]["COV"]["avg"] = aux2["media"]
            datos[d]["COV"]["med"] = aux2["mediana"]
        else :
            print("WARNING: {} no encontrado. Ejecuta AUP/coveragePanells.R para generarlo".format(cov))

        # Guardar en el diccionario el histograma con los cambios de base unica (SNVs) de cada muestra
        anno = "{}/variantCalling/raw.hg19_multianno.txt".format(path)
        kaks, fwrv, hist = vcfQC.getRatios(anno, False)
        datos[d]["VCF"] = hist

        # Guardar el path del archivo de coverage para el script de R que creara los graficos
        cov = "{}/coverage/coverageBase.txt".format(path)
        if os.path.isfile(cov) :
            datos[d]["COV_PATH"] = cov
        else :
            print("WARNING: Archivo de coverage ({}) no encontrado en la muestra {}".format(cov, d))

    return datos

def scriptR(datos) :
    filename = "globalQC.R"
    fqreads = "fq <- c("
    noms = "samps <- c("
    for d in datos.keys() :
        fqreads += "{},".format(datos[d]["FQ"])
        noms += "'{}',".format(d)

    fqreads = fqreads.rstrip(",")
    noms = noms.rstrip(",")
    fqreads += ")\n"
    noms += ")\n"

    with open(filename, "w") as fi :
        fi.write("library(RColorBrewer) # Paleta de colores para los graficos")
        fi.write("# Grafico de barras con la cantidad de reads de cada muestra\n")
        fi.write(fqreads)
        fi.write(noms)
        fi.write("png('FASTQ_distribution.png', width = 720, height = 720)\n")
        fi.write("barplot(fq, names.arg = samps, col = brewer.pal(12, 'Set3'))\n")
        fi.write("dev.off()")
        for k, v in datos.items() :
            fi.write("# Grafico de barras con el porcentaje de bases con un coverage determinado\n")
            fi.write("{smp} <- c()")



def main() :
    # Leer la carpeta en la que se ha invocado el programa para comprobar si se va a analizar la tanda en curso o hay que preguntar que tanda analizar
    dir = os.getcwd().split("/")[-1]
    dc = {}
    if dir.startswith(cte.prefijoTanda) :
        dc = recogerDatos()
    else :
        tanda = input("INPUT: Numero de tanda donde se realiza el control de calidad: ")
        try :
            tanda = int(tanda)
        except ValueError :
            print("ERROR: Numero de tanda erroneo")
            sys.exit(1)
        path = "{wd}/tanda{tn}".format(wd = cte.workindir, tn = tanda)
        dc = recogerDatos(path)

    print(dc)
    scriptR(dc)

    # Executar el Rscript en la cerpeta que toca

if __name__ == "__main__" :
    main()
