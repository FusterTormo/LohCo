#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math

#Convereix els bytes en una mesura mes llegible per humans
def convert_size(size_bytes):
    ret = ""
    if size_bytes == 0 :
        ret = "0 B"
    else :
        size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
        i = int(math.floor(math.log(size_bytes, 1024)))
        p = math.pow(1024, i)
        s = round(size_bytes / p, 2)
        ret = "{} {}".format(s, size_name[i])
    return ret

def eliminarVariantCalling(ruta) :
    print("INFO: Eliminant temporals de la carpeta de variants")
    cont = 0
    pes = 0
    for root, dirs, files in os.walk(ruta) :
        break

    for d in dirs :
        path = "{root}/{samp}/variantCalling".format(root = ruta, samp = d)
        arxius = os.listdir(path)
        for a in arxius :
            if a == "bwa.pileup" :
                pes += os.path.getsize("{dir}/bwa.pileup".format(dir = path))
                cont += 1
                os.remove("{dir}/bwa.pileup".format(dir = path))
            elif a.startswith("filtro") and not a.endswith(".allInfo") :
                pes += os.path.getsize("{dir}/{file}".format(dir = path, file = a))
                cont += 1
                os.remove("{dir}/{file}".format(dir = path, file = a))
    print("INFO: {} arxius eliminats. {} espai alliberat".format(cont, convert_size(pes)))

def eliminarAliniament(ruta) :
    print("INFO: Eliminant bams intermedis")
    cont = 0 # Contador de arxius elminats
    pes = 0
    arxius = ["bwa.sam", "bwa.sort.bam", "bwa.sort.bai", "indelsGATK.list"] # Lista d'arxius a eliminar
    # Recollir els directoris de cadascuna de les mostres en la tanda
    for root, dirs, files in os.walk(ruta) :
        break

    # En cadascuna de les mostres, eliminar els arxius bwa.sam, bwa.sort.bam, bwa.sort.bai i indelsGATK.list
    for d in dirs :
        path = "{root}/{samp}/bwaAlign".format(root = ruta, samp = d)
        for a in arxius :
            arx = "{}/{}".format(path, a)
            if os.path.isfile(arx) :
                pes += os.path.getsize(arx)
                os.remove(arx)
                cont += 1
    print("INFO: {} bams eliminats. {} espai alliberat".format(cont, convert_size(pes)))


def eliminarFASTQ(ruta) :
    print("INFO: Eliminant arxius FASTQ")
    cont = 0
    pes = 0
    for root, dirs, files in os.walk(ruta) :
        for f in files :
            if f.endswith("fastq.gz") :
                cont += 1
                pes += os.path.getsize("{path}/{arxiu}".format(path = ruta, arxiu = f))
                os.remove("{path}/{arxiu}".format(path = ruta, arxiu = f))

    print("INFO: {} arxius eliminats. {} espai alliberat".format(cont, convert_size(pes)))

def main(tanda = "") :
    root = "/home/ffuster/panalisi/resultats"
    if tanda == "" :
        tanda = input("REQUEST: Dona'm el numero de tanda que vols netejar: ")

    try :
        t = int(tanda)
        tanda = "tanda{}".format(t)
    except ValueError :
        print("ERROR: El numero de tanda no es un numero. Valor introduit: {}".format(tanda))
        sys.exit(1)

    path = "{root}/{tanda}".format(root = root, tanda = tanda)
    if os.path.isdir(path) :
        eliminarFASTQ(path)
        eliminarAliniament(path)
        eliminarVariantCalling(path)
    else :
        print("ERROR: La carpeta {} no existeix".format(path))

if __name__ == "__main__" :
    if len(sys.argv) > 1 :
        main(sys.argv[1])
    else :
        main()
