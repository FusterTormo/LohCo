#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Funciones para crear el log de analisis de los paneles
"""

"""
FUNCTIONS:
    prepararScript - Crea el bash con todos los comandos que se van a usar en el análisis del panel
    extraerRG - Extrae el Read Group de los FASTQ
    FALTEN FUNCIONS PER POSCAR ACI
"""
import os

"""
CONSTANTS:
    Rutas de los programas y parametros de cada uno de los pasos dentro de la pipeline
"""
fastqc = "/opt/FastQC/fastqc" #Ruta al FastQC (control de calidad de los FASTQ)
pfqc = "-o fastqc/ -f fastq -extract -q -t 6 {input}" # Parametros en la ejecucion dels fastqc
bwa = "bwa"
pbwa = ""
gatk = ""
picard = ""
vc = "" #Ruta al variant caller que se va a usar (Strelka2)
anno = "" #Ruta al ANNOVAR (anotador de variantes)
cov = "" #Script de coverage que se va a hacer

referencia = "/home/ffuster/panalisi/referencies/gatkHg19.fa"
manifest = "/home/ffuster/panalisi/resultats/manifest.bed"
indels = "/home/ffuster/panalisi/referencies/gold_indels.vcf"
dbsnp = "/home/ffuster/panalisi/referencies/dbsnp_138.hg19.vcf"
genes = "/home/ffuster/panalisi/resultats/gensAestudi.txt"

pathAnalisi = "/home/ffuster/panalisi/resultats" # Ruta donde se ejecutan y guardan los analisis
prefijoTanda = "tanda" # Prefijo que tiene todas las tandas analizadas

def extraerRG(fastq) :
    """
    Extrae el Read Group de un archivo FASTQ.

    Se espera el formato de genomica del IGTP: ID_SAMPLE_L001_R1.fastq.gz

    Parameters
    ----------
        fastq : str
            Nombre del archivo fastq del que se quiere sacar el read group
    Returns
    -------
        str
            La cadena de read group lista para incrustar en el comando de BWA.
    """
    pass

def getFASTQnames(path) :
    """
    Recoger los nombres de los archivos FASTQ que hay en la ruta especificada.

    Parameters
    ----------
        path : str
            Ruta absoluta donde estan los FASTQ que se van a analizar

    Returns
    -------
        list
            Lista con las rutas absolutas de los FASTQ encontrados en la ruta que se paso como parametro
    """
    files2copy = []
    print("INFO: Recollint el nom dels FASTQ des de {}".format(path))
    for root, dirs, files in os.walk(path) :
        for fic in files :
            # Parche. Se asume que el archivo FASTQ tiene extension .fastq.gz
            aux, extension = os.path.splitext(fic) # Eliminar la primera extension
            name, extension2 = os.path.splitext(aux)
            # Coger los FASTQ que se copiaran en la carpeta de la tanda
            if extension2 == ".fastq" :
                pt = "{}/{}".format(root, fic)
                files2copy.append(pt)
    print("INFO: {} arxius trobats".format(len(files2copy)))

    return files2copy

def getTanda() :
    """
    Busca en la carpeta de analisis cual es el numero de la ultima tanda. Devuelve el numero de la tanda siguiente

    Returns
    -------
        int
            Numero asignado para la siguiente tanda que se va a analizar
    """
    nums = []
    prefijo = len(prefijoTanda)
    for root, dirs, files in os.walk(pathAnalisi) :
        for d in dirs :
            if d.startswith(prefijoTanda) :
                aux = int(d[prefijo:])
                nums.append(aux)
        break
    sig = max(nums) + 1
    return sig

def doListaGenes() :
    """
    Crear el archivo con la lista de genes que hay dentro del manifest
    """
    listaGenes = []
    with open(manifest, "r") as fi :
        for l in fi :
            aux = l.split("\t")[3]
            aux = aux.strip()
            if aux not in listaGenes :
                listaGenes.append(aux)
    with open(genes, "w") as fi :
        fi.write("\n".join(listaGenes))

# TODO: Documentar correctament
def comprobarArchivos() :
    """
    Comprueba si los archivos necesarios para el analisis existen en la ruta especificada en las constantes
    """
    print("INFO: Buscant els arxius necessaris per executar la pipeline")
    if not os.path.isfile(referencia) :
        raise IOError("No se encuentra el genoma de referencia")
    if not os.path.isfile(manifest) :
        raise IOError("No se encuentra el manifest")
    if not os.path.isfile(indels) :
        raise IOError("No se encuentra el archivo para poder realizar el realineamiento de indels")
    if not os.path.isfile(dbsnp) :
        raise IOError("No se encuentra el archivo de SNPs")
    if not os.path.isfile(genes) :
        print("WARNING: No se encuentra el archivo con la lista de genes del manifest. Creando el archivo")
        doListaGenes()


def prepararScript(ruta) :
    """
    Programa principal de la libreria. Prepara el log con todos los comandos necesarios para lanzar la pipeline

    Comandos necesarios:
    * Averiguar el numero de tanda
    * Montar la estructura de la tanda
    * Copiar los FASTQ
    * Crear el bash de analisis
    * Crear, si no existe, la lista de los genes que contiene el manifest (gensAestudi.txt)
    * Crear el log con todos los comandos para cada una de las muestras del panel
    * Ejecutar, si procede el analisis usando la libreria subprocess

    Parameters
    ----------
        ruta : str
            Ruta absoluta donde esta la carpeta con los FASTQ que se van a analizar en esta pipeline
    """
    os.chdir(pathAnalisi) # Cambiar el directorio de trabajo a la carpeta de analisis
    tnd = getTanda() # Crear el nombre de la carpeta donde se guardaran los analisis y el nombre del bash con todos los comandos
    tanda = "{prefijo}{tanda}".format(prefijo = prefijoTanda, tanda = tnd)
    arxiu = "logTanda{tanda}.sh".format(tanda = tnd)

    print("INFO: Els resultats de l'analisi es guardaran en {path}/{tanda}".format(path = pathAnalisi, tanda = tanda))
    comprobarArchivos() # Esta funcion dispara una excepcion en caso de que no se encuentre alguno de los archivos necesarios para el analisis
    fastqs = getFASTQnames(ruta)
    if len(fastqs) == 0 :
        print("ERROR: No s'han trobat arxius FASTQ en {}".format(ruta))
        sys.exit(1)
    else :
        print("INFO: Creant el bash per la tanda {}".format(tnd))
        with open(arxiu, "w") as fi :
            fi.write("#!/bin/bash\n\n") # Shebang del bash
            fi.write("#Referencias usadas en este analisis\n")
            fi.write("ref={}\n".format(referencia))
            fi.write("mani={}\n".format(manifest))
            fi.write("indels={}\n".format(indels))
            fi.write("sites={}\n".format(dbsnp))
            fi.write("gens={}\n\n".format(genes))

            #Crear la funcion que copia los datos (FASTQs) en la carpeta de analisis
            fi.write("function copiar {\n")
            fi.write("\tcd {}\n".format(pathAnalisi))
            fi.write("\tmkdir {tanda} ; cd {tanda}\n\n".format(tanda = tanda))
            fi.write("\techo -e \"################################\\n\\tCopiant dades\\n################################\\n\"\n")
            for f in fastqs :
                patx = f.replace(" ", "\\ ") #Parche para leer los espacios en la terminal bash
                fi.write("\trsync -aP {} .\n".format(patx))

            fi.write("\tmv ../{} .\n".format(arxiu))
            fi.write("}\n\n")

            # TODO Crear les comandes per cadascuna de les etapes de l'analisi
            # TODO: Esta part es com si posara analizar.sh dins del log
            fi.write("function analisi {\n")
            fi.write("\tforward=$1\n\treverse=$2\n\treadgroup=$3\n\talias=$4\n")
            fi.write("\tmkdir $alias\n")
            fi.write("\tcd $alias\n")
            fi.write("\t$HOME/anpanmds/fastqc.sh ../$forward ../$reverse\n")
            fi.write("bamtobed per fer els analisis de ON TARGET")
            fi.write("\t$HOME/anpanmds/align.sh $reference ../$forward ../$reverse $readgroup\n")
            fi.write("\t$HOME/anpanmds/postAlign.sh recalibrate bwaAlign/bwa.sort.bam\n")
            fi.write("coverage")
            fi.write("estadistiques de l'analisi: on target, off target, % bases amb X coverage, resum dels tests, % duplicats (si cal), grafiques de coverage")
            fi.write("\t$HOME/anpanmds/variantCalling.sh $manifest bwaAlign/bwa.recalibrate.bam\n")
            fi.write("\T$HOME/anpanmds/variantAnnotation.sh variants.vcf")
            fi.write("Script per re-anotar")
            fi.write("Script per filtrar")
            fi.write("}\n\n")

        # TODO: Crear les comandes per analitzar de la mateixa manera a com s'esta fent en els panells d'ALL
