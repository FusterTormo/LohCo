#!/usr/bin/python
# -*- coding: utf-8 -*-

import subprocess
import sys
import os
import re

def pct(bam = "bwaAlign/bwa.nodup.bam") :
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
                with open(output = "alnQC.txt", "a") as fi :
                    fi.write("DUPS: {}".format(pct))
        else :
            print("ERROR: Samtools no se ejecuto correctamente. Descripcion: {}".print(err))
    else :
        print("ERROR: Bam no encontrado. Ruta de busqueda: {}".format(bam))



if __name__ == "__main__" :
    bam = sys.argv[1]
    pct(bam)
