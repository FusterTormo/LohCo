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
    # TODO: Eliminar este print
    print(dirs)
    print("INFO: {} muestras encontradas".format(len(dirs)))
    datos = {}
    for d in dirs :
        datos[d] = {"FQ" : 0, "BAM" : 0, "ON" : 0, "OFF"}
        path = "{}/{}".format(ruta, d)
        qc = "{}/{}".format(path, cte.qcaln)
        if os.path.isfile(qc)) :
            with open(qc, "r") as fi :
                aux = fi.read()
            aux2 = eval(aux)
            datos[d]["FQ"] = aux2["FASTQ"]
            datos[d]["BAM"] = aux2["BAM"]
            datos[d]["ON"] = aux2["ON"]
            datos[d]["OFF"] = aux2["OFF"]
        else :
            print("ERROR: {} no encontrado. Ejecuta AUP/bamQC.py para generarlo".format(qc))


def main() :
    # Leer la carpeta en la que se ha invocado el programa para comprobar si se va a analizar la tanda en curso o hay que preguntar que tanda analizar
    dir = os.getcwd().split("/")[-1]
    if dir.startswith(cte.prefijoTanda) :
        recogerDatos()
    else :
        tanda = input("INPUT: Numero de tanda donde se realiza el control de calidad: ")
        try :
            tanda = int(tanda)
        except ValueError :
            print("ERROR: Numero de tanda erroneo")
            sys.exit(1)
        path = "{wd}/tanda{tn}".format(wd = cte.workindir, tn = tanda)
        recogerDatos(path)
