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
        datos[d] = {"FQ" : 0, "BAM" : 0, "ON" : 0, "OFF" : 0, "COV" : {}, "VCF" : {}}
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
            datos[d]["COV"]["0"] = aux2["bases0"]
            datos[d]["COV"]["30"] = aux2["bases30"]
            datos[d]["COV"]["100"] = aux2["bases100"]
            datos[d]["COV"]["500"] = aux2["bases500"]
            datos[d]["COV"]["1000"] = aux2["bases1000"]
        else :
            print("WARNING: {} no encontrado. Ejecuta AUP/coveragePanells.R para generarlo".format(cov))

        # Guardar en el diccionario el histograma con los cambios de base unica (SNVs) de cada muestra
        anno = "{}/variantCalling/raw.hg19_multianno.txt".format(path)
        kaks, fwrv, hist = vcfQC.getRatios(anno, False)
        datos[d]["VCF"] = hist

    return datos


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

if __name__ == "__main__" :
    main()
