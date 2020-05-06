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
    """
    Convierte el manifest en un diccionario

    Abre el archivo manifest y convierte todas las regiones en un diccionario con el formato: {"chr" : [[ini_reg1 - fin_reg1], [ini_reg2 - fin_reg2]]}

    Returns
    -------
        dict
            El manifest convertido en diccionario
    """
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
    """
    Comprueba si un read pasado por parametro esta dentro del manifest o no

    Se comprueba si, al menos, una parte del read esta dentro de las regiones del manifest.

    Parameters
    ----------
        read : str
            Cadena de texto con el formato cromosoma\tinicio\tfin
        manifest : dict
            Diccioario en el que buscar si el read esta dentro o no

    Returns
    -------
        bool
            True si el read esta dentro del manifest. False en caso contrario
    """
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
    """
    Programa principal

    Recoge los datos de calidad del bam. Se consideran parametros de calidad el numero de reads que se tenian al principio. El numero de reads que contiene el bam. El numero de reads que estan dentro
    del manifest (ON target) y el numero de reads fuera del manifest (OFF target). Un read se considera ON target si tiene, al menos, una base dentro de las regiones definidas por el manifest.
    Se guardan todos los datos de calidad en un archivo de texto llamado alnQC.txt
    """
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
            print("ERROR: Bed no encontrado en ruta {}".format(bed))
    else :
        print("ERROR: Manifest no encontrado en ruta {}".format(manifest))

    fqreads = sum(reads)
    breads = on + off
    print("INFO: Escribiendo resultados en {}".format(output))
    with open(output, "w") as fi :
        fi.write("{")
        fi.write("\'FASTQ\': {},".format(fqreads))
        fi.write("\'BAM': {},".format(breads))
        fi.write("\'ON\': {},".format(on))
        fi.write("\'OFF\': {}".format(off))
        fi.write("}")

if __name__ == "__main__" :
    main()
