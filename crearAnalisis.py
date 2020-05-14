# -*- coding: utf-8 -*-
#!/usr/bin/python

import os
import re
import subprocess
import sys

#Librerias propias
import getCommands as gc

# Constantes
hgref = "/home/ffuster/share/biodata/solelab/referencies/gatk_hg38/hg38.fa"
snpSites = "/home/ffuster/share/biodata/solelab/referencies/gnomad.genomes.r3.0.sites.vcf"
logfile = "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/analisi_Fco/log.sh"
wd = "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/analisi_Fco"

def crearRG(fq) :
    """Crear el read group para el FASTQ pasado por parametro"""
    aux = fq.split("_")
    # El formato (arbitrario) del read group sera ID=flowcell, sample=index
    rg = "\"@RG\\tID:{id}\\tSM:S{samp}\\tPL:ILLUMINA\"".format(id = aux[0], samp = aux[1])
    return rg

def leerMuestras() :
    """
    Lee los excels que contienen la informacion de las muestras. Separar las muestras con formato: {"id" : {"fastq" , [[fw1, rv1, rg1], [fw2, rv2, rg2], ...}, "gender" : "[M]/[F]"}
    El diccionario de muestras control no tendra informacion del genero. Pero se tiene en cuenta que, de los controles generales, C2 es una mujer y C3 es un hombre
    """
    arx = ["MDS_12.csv", "MDS_13.csv", "MDS_14.csv", "MDS_15.csv"]
    paths = {"MDS_12.csv" : "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/MDS_12/20190618/FASTQ", "MDS_13.csv" : "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/MDS_13/20190612/FASTQ",
    "MDS_14.csv" : "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/MDS_14/20200114/FASTQ", "MDS_15.csv" : "/home/ffuster/share/gsole/WES_Projecte_SMD_NMP/MDS_15/20200114/FASTQ"}
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
                    rg = crearRG(fq1.split("/")[-1])
                    if id.startswith("CD3") or id.endswith("CD3") : # El identificador indica que es un control
                        if id not in controls.keys() :
                            controls[id] = {"fastq" : [[fq1, fq2, rg]]}
                        else :
                            controls[id]["fastq"].append([fq1, fq2, rg])
                    else :
                        if id not in tumors.keys() :
                            tumors[id] = {"fastq" : [[fq1, fq2, rg]], "gender" : aux[10]}
                        else :
                            tumors[id]["fastq"].append([fq1, fq2, rg])
    return tumors, controls

def getControl(muestra, controles) :
    # Extraer el numero identificador de la muestra tumoral
    regex = "[\d]+"
    match = re.search(regex, muestra)
    idmuestra = ""
    if match :
        id = match.group()
        for c in controles.keys() :
            if c.endswith(id) :
                idmuestra = c
                break
    else :
        print("WARNING: Identificador no encontrado en {}".format(muestra))

    return idmuestra

def ejecutar(orden) :
    proc = subprocess.Popen(orden, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = proc.communicate()
    if proc.returncode == 0:
        # Escribir la orden finalizada en el log.
        with open(logfile, "a") as fi :
            fi.write("{}\n".format(orden))
    else :
        print("ERROR: Al ejecutar {cmd}\n\nDescripcion: {dsc}\nComando: {cmd}".format(cmd = orden, dsc = err))
        sys.exit()

def alinear(muestra) :
    """Alinea, ordena, recalibra y marca duplicados para la muestra pasada por parametro"""
    sams = []
    # Alinear
    if len(muestra["fastq"]) > 1 :
        it = 1
        for sample in muestra["fastq"] :
            aux = "bwa_{}.sam".format(it)
            sams.append(aux)
            cmd = gc.getAln(sample[2], hgref, sample[0], sample[1], aux)
            ejecutar(cmd)
            it += 1
    else :
        sams.append("bwa.sam")
        cmd = gc.getAln(sample[2], hgref, sample[0], sample[1], "bwa.sam")
        ejecutar(cmd)
    # Ordenar
    if len(sams) == 1 :
        cmd = gc.getPcSort("bwa.sam", "bwa.sort.bam")
    else :
        cmd = gc.getPcMerge(sams, "bwa.sort.bam")
    ejecutar(cmd)
    # A partir de este paso ya se trabaja con un unico archivo bam
    ejecutar(gc.getPcIndex("bwa.sort.bam"))
    ejecutar(gc.getGATKrecal("bwa.sort.bam", hgref, snpSites, "bwa.recal.bam"))
    ejecutar(gc.getPcMarkduplicates("bwa.recal.bam", "bwa.nodup.bam"))
    # Eliminar todos los archivos temporales
    # for a in sams :
    #     os.remove(a)
    # os.remove("bwa.sort.bam")
    # os.remove("bwa.sort.bai")
    # os.remove("bwa.recal.bam")
    # os.remove("bwa.recal.bai")

def main() :
    tm, cn = leerMuestras()
    for k,v in tm.items() :
        c = getControl(k, cn)
        if c == "" :
            if v["gender"] == "F" :
                c = "C2"
            elif v["gender"] == "M" :
                c = "C3"
        if c != "" and c != "C2" and c != "C3" :
            folder = "{}vs{}".format(k,c)
            if os.path.isdir(folder) : # La carpeta existe, puede que algo se haya ejecutado
                #Guardar todos los archivos de la carpeta en una lista
                filenames = os.listdir(folder)
                os.chdir(folder)
                if "bwa.nodup.bam" not in filenames :
                    alinear(v)
                else :
                    print("TODO: Comprobar si los variant calls se han hecho en la carpeta")
            else :
                print("INFO: Analizando {}".format(folder))
                os.makedirs(folder, 0o754)
                # Ordenes para alinear la muestra
                os.chdir(folder)
                with open(logfile, "a") as fi :
                    fi.write("cd {}\n".format(os.path.abspath(os.getcwd())))
                alinear(v)
                os.chdir(wd)
                sys.exit()
            #alinear(tm)
            #alinear(cn)

if __name__ == "__main__" :
    main()
