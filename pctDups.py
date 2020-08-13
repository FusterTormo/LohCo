#!/usr/bin/python
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import re

output = "alnQC.txt"

def pct(bam = "bwaAlign/bwa.nodup.bam") :
    """
    Extraer el porcentaje de duplicados de un bam pasado por parametro

    Ejecuta samtools flagstat para extraer el numero de reads considerados duplicados y luego calcula el porcentaje de duplicados que contiene el bam pasado por parametro.
    El resultado se guarda en un archivo de texto.

    Parameters
    ----------
        bam : str, optional
            Ruta donde se encuentra el bam en el cual se quiere calcular el porcentaje de duplicados
    """
    cmd = "samtools flagstat {}".format(bam)
    dups = -1
    reads = -1
    regex = " in total "
    if os.path.isfile(bam) :
        pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        std, err = pr.communicate()
        if proc.returncode == 0 :
            for i in std.split("\n") :
                if i.endswith("duplicates") :
                    dups = float(i.split(" + ")[0])
                if re.search(regex, i) :
                    reads = float(i.split(" + ")[0])

            if reads > -1 and dups > -1 :
                pct = dups/reads
                pct *= 100
                with open(output, "r") as fi :
                    aux = fi.read().strip("}")
                aux += ",\'DUPS\': {}".format(pct)
                aux += "}"
                with open(output, "w") as fi :
                    fi.write(aux)
        else :
            print("ERROR: Samtools no se ejecuto correctamente. Descripcion: {}".format(err))
    else :
        print("ERROR: Bam no encontrado. Ruta de busqueda: {}".format(bam))

if __name__ == "__main__" :
    bam = sys.argv[1]
    pct(bam)
