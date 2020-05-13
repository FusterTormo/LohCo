# -*- coding: utf-8 -*-
#!/usr/bin/python

import pandas

#Librerias propias
import getCommands

def leerMuestras(ruta) :
    """
    Lee un archivo excel pasado por parametro, separa las muestas y crea la estructura
    { "idMuestra" : {"tumor" : [ fw1.fastq, rv1.fastq, fw2.fastq, rv2.fastq ...], "control" : [ fw1.fastq, rv1.fastq ]}}
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
                            tumors[id] = {"fastq" : [[fq1, fq2]]}
                        else :
                            tumors[id]["fastq"].append([fq1, fq2])
