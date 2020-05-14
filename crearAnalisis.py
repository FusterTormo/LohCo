# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
import subprocess

#Librerias propias
import getCommands

def leerMuestras() :
    """
    Lee los excels que contienen la informacion de las muestras. Separar las muestras con formato: {"id" : {"fastq" , [[fw1, rv1, rg1], [fw2, rv2, rg2], ...}, "gender" : "[M]/[F]"}
    El diccionario de muestras control no tendra informacion del genero. Pero se tiene en cuenta que, de los controles generales, C2 es una mujer y C3 es un hombre
    """
    arx = ["MDS_12.csv", "MDS_13.csv", "MDS_14.csv", "MDS_15.csv"]
    paths = {"MDS_12.csv" : "../MDS_12/20190618/FASTQ", "MDS_13.csv" : "../MDS_13/20190612/FASTQ", "MDS_14.csv" : "../MDS_14/20200114/FASTQ", "MDS_15.csv" : "../MDS_15/20200114/FASTQ"}
    tumors = {}
    controls = {}
    for a in arx :
        cabecera = True
        with open(a, "r") as fi :
            dirpath = paths[a]
            for l in fi :
                if cabecera :
                    cabecera = False
                else :
                    aux = l.strip("\n").split(";")
                    id = aux[6]
                    fcell = aux[0].split("(")[0]
                    fq1 = "{dir}/{flowcell}_{lane}_{index}_1.fastq.gz".format(dir = dirpath, flowcell = fcell, lane = aux[1], index = aux[2])
                    fq2 = "{dir}/{flowcell}_{lane}_{index}_2.fastq.gz".format(dir = dirpath, flowcell = fcell, lane = aux[1], index = aux[2])
                    if id.startswith("CD3") or id.endswith("CD3") : # El identificador indica que es un control
                        if id not in controls.keys() :
                            controls[id] = {"fastq" : [[fq1, fq2]]}
                        else :
                            controls[id]["fastq"].append([fq1, fq2])
                    else :
                        if id not in tumors.keys() :
                            tumors[id] = {"fastq" : [[fq1, fq2]], "gender" : aux[10]}
                        else :
                            tumors[id]["fastq"].append([fq1, fq2])
    return tumors, controls

def getControl(muestra, controles) :
    # Extraer el numero identificador de la muestra tumoral
    regex = "[\d]+"
    match = re.search(regex, muestra)
    idmuestra = ""
    if match :
        id = match.group()
        print("INFO: Encontrado {} en la muestra {}".format(id, muestra))
        for c in controles.keys() :
            if c.endswith(id) :
                idmuestra = c
                break
    else :
        print("WARNING: Identificador no encontrado en {}".format(muestra))

    return idmuestra


def main() :
    tm, cn = leerMuestras()
    for k,v in tm.items() :
        c = getControl(k, cn)
        folder = "{}_VS_{}".format(t,c)
        print(folder)
    #alinear(tm)
    #alinear(cn)

if __name__ == "__main__" :
    leerMuestras()
