#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Comprobar la calidad del alineamiento
"""

import os

# IDEA: Algunes d'aquestes constants podrien anar en una llibreria de constants
bed = "bwaAlign/bwa.bed"
fastqc = "fastqc"
fastqcFI = "fastqc_data.txt"
manifest = "/home/ffuster/panalisi/resultats/manifest.bed"
output = "alnQC.txt"

def convertirManifest() :
    m = {}
    with open(manifest, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[0]
            ini = int(aux[1])
            fin = int(aux[2])
            if chr in m.keys() :
                m[chr].append([ini, fin])
            else :
                m[chr] = [[ini, fin]]
    return m

def enManifest(read, manifest) :
    esta = False
    aux = read.split("\t")
    chr = aux[0]
    ini = int(aux[1])
    fin = int(aux[2])
    if chr in manifest.keys() :
        for r in manifest[chr] :
            if fin >= r[0] and ini <= r[1] :
                esta = True
                break
    return esta

def main() :
    on = 0
    off = 0
    fqreads = 0
    bamreads = 0
    reads = []

    # Comprobar si existe la carpeta de FASTQC. Extraer la informacion necesaria de los FASTQ en tal caso
    print("INFO: Leyendo FASTQ")
    if os.path.isdir(fastqc) :
        for dir, dirnames, filenames in os.walk(fastqc) :
            break

        for d in dirnames :
            with open("{}/{}/{}".format(fastqc, d, fastqcFI), "r")  as fi :
                for l in fi :
                    if l.startswith("Total Sequences") :
                        reads.append(int(l.split("\t")[-1].strip()))
                        break

    # Comprobar si existe el manifest. Convertirlo en un diccionario para acceder a los reads mas rapidamente
    print("INFO: Leyendo manifest")
    if os.path.isfile(manifest) :
        dic = convertirManifest()
        print("INFO: Leyendo alineamiento")
        if os.path.isfile(bed) :
            with open(bed, "r") as fi :
                for l in fi :
                    if enManifest(l, dic) :
                        on += 1
                    else :
                        off += 1
        else :
            print("ERROR: Bed no encontrado. Buscando en ruta {}")
    else :
        print("ERROR: Manifest no encontrado. Buscando en ruta {}")

    fqreads = sum(reads)
    breads = on + off
    print("INFO: Escribiendo resultados en {}".format(output))
    with open(output, "w") as fi :
        fi.write("FASTQ: {}\n".format(fqreads))
        fi.write("BAM: {}\n".format(breads))
        fi.write("ON: {}\n".format(on))
        fi.write("OFF: {}\n".format(off))

if __name__ == "__main__" :
    main()
