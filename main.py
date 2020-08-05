#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Programa principal. Muestra los pasos a seguir para poder hacer el analisis automatico de paneles en SMD. Usa la libreria creaLog para ejecutar el analisi
"""

import os
import subprocess

import bamQC as bq
import creaLog as cl
import data2excel as xls
import filtStrelka as fis
import filtMutect as fim
import getCommands as gc
import manifestOp as op

def fastq() :
    """Analizar una unica muestra"""
    # TODO: Recoger la ruta absoluta de los FASTQ. Lanzar custom para ver el tipo de analisis a ejecutar
    pass

def reanalizar() :
    """Reanalizar una muestra sin borrar lo que ya esta guardado"""
    # # TODO: COMO LO HAGO ??????????
    pass

def vcf() :
    """Anotar y filtrar un vcf. Guardar los datos en un excel

    Pide al usuario la ruta de un archivo vcf. Anota y filtra los datos de dicho excel. Finalmente guarda todo en un excel. Todos los datos se guardan en la carpeta donde esta el vcf
    """
    ruta = input("Introducir el path absoluto del vcf: ")
    name = input("Introducir el nombre de la muestra: ")
    hg = input("Introducir el genoma de referencia usado (hg19, hg38): ")
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
                    vcal = "Sterlka2"
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
        ruta = input("Introducir el path absoluto del bam a analizar: ")

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
                    if input("Convertir los datos en excel? (S/N)").lower() == "s" :
                        ordenes.append("excel")
    print("INFO: Se van a ejecutar los siguientes pasos: {}".format(", ".join(ordenes)))
    ruta = input("Introducir el path absoluto de la carpeta donde estan los FASTQ a analizar: ")
    cl.prepararPanel(ruta, ordenes)

def anotarManifest() :
    ruta = input("Introducir ruta absolute del manifest a anotar: ")
    op.anotarManifest(ruta)

def GUI() :
    """
    Menu de interaccion con el usuario. Muestra las opciones de analisis disponibles
    """
    # TODO: Pintar-ho bonico
    # IDEA: Crear un menu amb els passos que es volen seguir per fer l'analisi mes personalitzat
    ruta = input("Introducir el path absoluto de la carpeta donde estan los FASTQ a analizar: ")
    cl.prepararPanel(ruta, ["copiar","fastqc", "aln", "recal", "bamqc", "coverage", "mutectGerm", "vanno", "filtrar", "excel"])
    # TODO: Preguntar si executar la pipeline usant subprocess


if __name__ == "__main__" :
    GUI()
