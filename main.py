#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
MAIN: Programa principal. Muestra los pasos a seguir para poder hacer el analisis automatico de paneles en SMD. Usa la libreria creaLog para ejecutar el analisi
"""

import os
import subprocess
import sys

import bamQC as bq
import creaLog as cl
import data2excel as xls
import filtStrelka as fis
import filtMutect as fim
import getCommands as gc
import manifestOp as op

def reanalizar() :
    """Reanalizar una muestra sin borrar lo que ya esta guardado"""
    workindir = "/home/ffuster/panalisi/resultats"
    samp = input("INPUT: Identificador de la muestra a reanalizar: ")
    # Buscar la muestra en el directorio de trabajo usando el comando find de bash
    cmd = "find {wd} -name {smp}*fastq.gz".format(smp = samp, wd = workindir)
    proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = proc.communicate()
    res = out.decode()
    if res == "" :
        print("INFO: No se ha encontrado la muestra en {}".format(workindir))
    else :
        # Convertir los resultados a una lista
        aux = res.split("\n")
        # Si se encuentran, recoger el directorio donde esta el fastq y lanzar un analisis con dicha muestra
        ruta = os.path.dirname(aux[0])
        acciones = custom()
        lanzarPanel(ruta, acciones)
        # Buscar el log que se ha creado al lanzar el panel
        for root, dirs, filenames in os.walk(workindir, topdown = True, ) :
            break

        for f in filenames :
            if f.startswith("logTanda") :
                logf = f
                break
        # Modificar el log para eliminar las otras muestras que podria haber en el directorio
        os.chdir(workindir)
        with open(logf, "r") as fi :
            cont = fi.readlines()
        for c in range(len(cont)) :
            if cont[c].startswith("mkdir") :
                cont[c] = "mkdir {smp} ; cd {smp}".format(smp = samp)
            if cont[c].startswith("mv") :
                cont[c] = "mv ../log{}.sh .".format(samp)
                print(cont[c])

        # Guardar los cambios en el archivo
        with open(logf, "w") as fi :
            fi.write(cont)
        os.rename(logf, "log{}.sh".format(samp)) # Cambiar el nombre logTanda__.sh por el nombre de la muestra

        print("WARNING: Si el directorio {} tiene los fastq de mas de una muestra, se habran creado analisis para todas las muestras. Elimina los analisis a las muestras que no interesan".format(ruta))

def vcf() :
    """Anotar y filtrar un vcf. Guardar los datos en un excel

    Pide al usuario la ruta de un archivo vcf. Anota y filtra los datos de dicho excel. Finalmente guarda todo en un excel. Todos los datos se guardan en la carpeta donde esta el vcf
    """
    ruta = input("INPUT: Introducir el path absoluto del vcf: ")
    name = input("INPUT: Introducir el nombre de la muestra: ")
    hg = input("INPUT: Introducir el genoma de referencia usado (hg19, hg38): ")
    # Comprobar que tipo de variant calling se ha hecho
    vcal = None
    format = "vcf4"
    # Buscar Mutect2, VarScan2, configureStrelka para saber cual de los variant callers se ha introducido
    with open(ruta, "r") as fi :
        for l in fi :
            # Leer solo la cabecera del vcf para determinar el variant caller y si el vcf es de comparacion somatica o somatico/germinal sin comparacion
            if l.startswith("#") :
                # Comprobar que la muestra es somatica y hace falta cambiar el parametro -format de ANNOVAR a vcf4old
                if l.startswith("#CHROM") :
                    if len(l.split("\t")) == 11 :
                        format = "vcf4old"

                if "Mutect2" in l :
                    vcal = "Mutect2"
                elif "VarScan2" in l :
                    vcal = "VarScan2"
                elif "configureStrelka" in l :
                    vcal = "Strelka2"
            else :
                break
    if vcal == None :
        print("ERROR: No se pudo determinar el variant caller del vcf")
    else :
        print("INFO: Anotando VCF")
        lst = None
        # Cambiar el directorio de trabajo a la carpeta donde esta el vcf
        wd = os.path.dirname(ruta)
        os.chdir(wd)
        # Anotar el vcf. Comprobar que no se haya anotado antes
        arx = "{}.{}_multianno.txt".format(name, hg)
        if not os.path.isfile(arx) :
            if hg == "hg19" :
                lst = gc.getANNOVAR(ruta, name, "hg19", format)
            elif hg == "hg38" :
                lst = gc.getANNOVAR(ruta, name, "hg38", format)
            else :
                raise ValueError("ERROR: Genoma de referencia no soportado")

            if lst != None :
                for c in lst.split("\n") :
                    proc = subprocess.Popen(c, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                    out, err = proc.communicate()

        # Ejecutar filtMutect, filtStrelka, dependiendo del variant caller encontrado
        if vcal == "Strelka2" :
            fis.main(arx, name)
            xls.crearExcel(name)
        elif vcal == "Mutect2" :
            fim.main(arx, name)
            xls.crearExcel(name)
        elif vcal == "VarScan2" :
            print("Pendiente hacer filtrado de VarScan2")

def bamQC(ruta = None) :
    """Calcula los parametros de calidad de un bam"""
    if ruta == None :
        ruta = input("INPUT: Introducir el path absoluto del bam a analizar: ")

    os.chdir(os.path.dirname(ruta))
    bq.bed = "bwa.bed"
    bq.bam = "bwa.recal.bam"
    bq.fastqc = "../fastqc"
    bq.main()
    # El script bamQC.py crea un archivo de texto, llamado alnQC.txt, con un JSON con los datos de interes
    with open("alnQC.txt", "r") as fi :
        txt = fi.read()
    qc = eval(txt)
    for k, v in qc.items() :
        if k == "ON" or k == "OFF":
            print("{} - {} ({:.2f} %)".format(k, v, 100*float(v)/float(qc["BAM"])))
        else :
            print("{} - {}".format(k, v))
    os.remove("bwa.bed")
    os.remove("alnQC.txt")

def custom() :
    """Crear un analisis customizado usando un wizard"""
    """
    OPCIONES YA CONTEMPLADAS EN CREAR LOG
    copiar
    fastqc
    aln
    recal
    bamqc
    mdups
    coverage
    strelkagerm
    mutectgerm
    vanno
    filtrar
    excel
    """
    ordenes = []
    print("\n\t------------------------------------------------\n\t\tAnalisis personalizado del panel\n\t------------------------------------------------\n")

    if input("Copiar los fastq? (S/N) ").lower() == "s" :
        ordenes.append("copiar")

    if input("Control de calidad de los fastq? (S/N) ").lower() == "s" :
        ordenes.append("fastqc")

    if input("Alinear? (S/N) ").lower() == "s" :
        ordenes.append("aln")
        if input("Recalibrar bases? (S/N) ").lower() == "s" :
            ordenes.append("recal")
        if input("Marcar duplicados? (S/N) ").lower() == "s" :
            ordenes.append("mdups")
        if input("Control de calidad del bam? (S/N) ").lower() == "s" :
            ordenes.append("bamqc")
        if input("Calculo de coverage del bam? (S/N) ").lower() == "s" :
            ordenes.append("coverage")
        aux = input("Con que programa hacer el variant calling? Opciones disponibles (strelka_germinal, mutect_germinal). En caso de no querer variant calling, escribir 'n' ")
        if aux.lower() != 'n' :
            aux = aux.replace("_germinal", "germ")
            ordenes.append(aux)
            if input("Anotar variantes usando ANNOVAR? (S/N) ").lower() == "s" :
                ordenes.append("vanno")
                if input("Filtrar variantes? (S/N) ").lower() == "s" :
                    ordenes.append("filtrar")
                    if input("Convertir los datos en excel? (S/N) ").lower() == "s" :
                        ordenes.append("excel")
    print("INFO: Se va a crear un log que ejecutara los siguientes pasos: {}".format(", ".join(ordenes)))
    return ordenes

def anotarManifest() :
    ruta = input("Introducir ruta absolute del manifest a anotar: ")
    op.anotarManifest(ruta)

def lanzarPanel(ruta, opciones) :
    cl.prepararPanel(ruta, opciones)

def GUI() :
    """
    Menu de interaccion con el usuario. Muestra las opciones de analisis disponibles
    """
    #system.clear()
    print("\t\t------------------------------------------------------\n\t\t\tAnalisis aUtomatico de Paneles\n\t\t------------------------------------------------------\n\nOpciones")
    print("1. Analisis tipico del panel")
    print("2. Analisis custom del panel")
    print("3. Analizar una unica muestra")
    print("4. Control de calidad de un bam")
    print("5. Anotar y filtrar un vcf")
    print("6. Reanalizar una muestra de una tanda, conservando el analisis previo")
    print("7. Anotar un manifest\n")
    opt = input("INPUT: Numero de opcion: ")
    # Comprobar si la opcion es un numero
    try :
        opt = int(opt)
    except ValueError :
        print("ERROR: La opcion no es un numero")
        sys.exit(1)
    # Analizar un panel de forma clasica/estandar
    if opt == 1 :
        ruta = input("INPUT: Introducir el path absoluto de la carpeta donde estan los FASTQ a analizar: ")
        lanzarPanel(ruta, ["copiar","fastqc", "aln", "recal", "bamqc", "coverage", "mutectGerm", "vanno", "filtrar", "excel"])
    # Analizar el panel usando una pipeline escogida por el usuario
    elif opt == 2 :
        ruta = input("INPUT: Introducir el path absoluto de la carpeta donde estan los FASTQ a analizar: ")
        opciones = custom()
        lanzarPanel(ruta, opciones)
    elif opt == 3 :
        ruta = input("INPUT: Introducir la carpeta donde estan los fastq de la muestra a analizar: ")
        opciones = custom()
        lanzarPanel(ruta, opciones)
    elif opt == 4 :
        ruta = input("INPUT: Introducir el path absoluto del bam a analizar: ")
        bamQC(ruta)
    elif opt == 5 :
        vcf()
    elif opt == 6 :
        reanalizar()
    elif opt == 7 :
        ruta = input("INPUT: Introducir la ruta del manifest: ")
        op.anotarManifest(ruta)
    else :
        print("ERROR: Opcion no valida")
        sys.exit(1)

if __name__ == "__main__" :
    GUI()
