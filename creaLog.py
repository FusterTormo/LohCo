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
import re
import sys
import getCommands as cmd

"""
CONSTANTS:
    Rutas de los programas y parametros de cada uno de los pasos dentro de la pipeline
"""
referencia = "/home/ffuster/share/biodata/solelab/referencies/ucsc/hg19.fa"
manifest = "/home/ffuster/panalisi/resultats/manifest.bed"
gzmanifest = "/home/ffuster/panalisi/resultats/manifest.bed.gz"
manifestidx = "/home/ffuster/panalisi/resultats/manifest.bed.gz.tbi"
"""Comandos para generar los manifest tal y como los necesitan los programas
    mv manifest.bed manifest_original.bed # Mantener una copia del original
    sort -k1,1 -k2,2n manifest_original.bed > manifestaux.bed
    sed 's/chr//g' manifestaux.bed > manifest.bed # Eliminar los 'chr'
    bgzip -c manifest.bed > manifest.bed.gz # Comprimir el archivo para Strelka2. Puede que la ruta original de bgzip sea /opt/htslib-1.10.2/bin/bgzip
    tabix -p bed manifest.bed.gz # Crear un indice para el manifest. Puede que la ruta original de tabix sea /opt/htslib-1.10.2/bin/tabix
"""
# Descargado desde https://gnomad.broadinstitute.org/downloads
dbsnp = "/home/ffuster/share/biodata/solelab/referencies/gnomad.exomes.r2.1.1.sites.vcf"
genes = "/home/ffuster/panalisi/resultats/gensAestudi.txt"

pathAnalisi = "/home/ffuster/panalisi/resultats" # Ruta donde se ejecutan y guardan los analisis
prefijoTanda = "tanda" # Prefijo que tiene todas las tandas analizadas
wd = "~/AUP" # Directorio de trabajo donde estan los scripts

def extraerRG(fastq) :
    """
    Extrae el Read Group de un archivo FASTQ.

    Extrae el identificador de muestra y el numero de muestra a partir del nombre de un archivo FASTQ pasado por parametro. Con estos datos construye los parametros necesarios
    para lanzar el analisis de la muestra. Estos parametros son: nombre del FASTQ R1, nombre del FASTQ R2, Read Group, Alias para la muestra
    Se espera el formato de genomica del IGTP: ID_SAMPLE_L001_R1.fastq.gz

    Parameters
    ----------
        fastq : str
            Nombre del archivo fastq del que se quiere sacar el read group
    Returns
    -------
        str
            La cadena de read group lista para incrustar en el comando de BWA. Si no se reconoce el formato, devolvera "No valido"
        str
            El identificador de la muestra. Si no se reconoce el formato, devolvera "No valido"
    """
    # Formato esperado ID_SAMPLE_L001_R1.fastq.gz
    aux = fastq.split("_")
    reg = "_S[\d]+_L001_R[\d]?_001\.fastq\.gz"
    if re.search(reg, fastq) :
        rgid = aux[0]
        sample = aux[1]
        fq1 = "{id}_{samp}_L001_R1_001.fastq.gz".format(id = rgid, samp = sample)
        fq2 = "{id}_{samp}_L001_R2_001.fastq.gz".format(id = rgid, samp = sample)
        rg = "\"@RG\\tID:{id}\\tSM:{samp}\\tPL:ILLUMINA\"".format(id = rgid, samp = sample)
        whole = "{fq1} {fq2} {rg} {alias}".format(fq1 = fq1, fq2 = fq2, rg = rg, alias = rgid)
    else :
        whole = "No valido"
        rgid = "No valido"
    return whole, rgid

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
    print("INFO: Recogiendo el nombre de los FASTQ en {}".format(path))
    for root, dirs, files in os.walk(path) :
        for fic in files :
            # Parche. Se asume que el archivo FASTQ tiene extension .fastq.gz
            aux, extension = os.path.splitext(fic) # Eliminar la primera extension
            name, extension2 = os.path.splitext(aux)
            # Coger los FASTQ que se copiaran en la carpeta de la tanda
            if extension2 == ".fastq" :
                pt = "{}/{}".format(root, fic)
                files2copy.append(pt)
    print("INFO: {} archivos encontrados".format(len(files2copy)))

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
    print("INFO: Buscando los archivos necesarios para ejecutar la pipeline")
    if not os.path.isfile(referencia) :
        raise IOError("No se encuentra el genoma de referencia")
    if not os.path.isfile(manifest) :
        raise IOError("No se encuentra el manifest")
    else : #Strelka2 quiere que el bam este comprimido e indexado. Comprobar si estas archivos existen
        if not os.path.isfile(gzmanifest) :
            raise IOError("ERROR: No se encuentra el bgzip del manifest, necesario para Strelka2. Ejecuta manifest.sh para crearlo")
        if not os.path.isfile(manifestidx) :
            raise IOError("ERROR: No se encuentra el indice del manifest, necesario para Strelka2. Ejecuta manifest.sh para crearlo")
    if not os.path.isfile(dbsnp) :
        raise IOError("No se encuentra el archivo de SNPs")
    if not os.path.isfile(genes) :
        print("WARNING: No se encuentra el archivo con la lista de genes del manifest. Creando el archivo")
        doListaGenes()

def prepararPanel(ruta, acciones) :
    """
    Programa principal de la libreria. Prepara el log con todos los comandos necesarios para lanzar la pipeline

    Comandos necesarios:
    * Averiguar el numero de tanda
    * Montar la estructura de la tanda
    * Copiar los FASTQ
    * Crear el bash de analisis
    * Crear, si no existe, la lista de los genes que contiene el manifest (gensAestudi.txt)
    * Crear el log con todos los comandos para cada una de las muestras del panel

    Parameters
    ----------
        ruta : str
            Ruta absoluta donde esta la carpeta con los FASTQ que se van a analizar en esta pipeline
        acciones : list
            Lista de acciones que se pueden agregar al bash. Valores aceptados: ["copiar","fastqc", "aln", "recal", "mdups", "bamqc", "coverage", "strelkaGerm", "mutectGerm", "vanno", "filtrar", "excel"]
    """
    os.chdir(pathAnalisi) # Cambiar el directorio de trabajo a la carpeta de analisis
    tnd = getTanda() # Crear el nombre de la carpeta donde se guardaran los analisis y el nombre del bash con todos los comandos
    tanda = "{prefijo}{tanda}".format(prefijo = prefijoTanda, tanda = tnd)
    arxiu = "logTanda{tanda}.sh".format(tanda = tnd)

    print("INFO: Los resultados del analisis se guardaran en {path}/{tanda}".format(path = pathAnalisi, tanda = tanda))
    comprobarArchivos() # Esta funcion dispara una excepcion en caso de que no se encuentre alguno de los archivos necesarios para el analisis
    fastqs = getFASTQnames(ruta)
    if len(fastqs) == 0 :
        print("ERROR: No se han encontrado archivos FASTQ en {}".format(ruta))
        sys.exit(1)
    else :
        print("INFO: Creando el bash para la tanda {}".format(tnd))
        with open(arxiu, "w") as fi :
            fi.write("#!/bin/bash\n\n") # Shebang del bash
            fi.write("#Referencias usadas en este analisis\n")
            fi.write("ref={}\n".format(referencia))
            fi.write("mani={}\n".format(manifest))
            fi.write("gzmani={}\n".format(gzmanifest))
            fi.write("sites={}\n".format(dbsnp))
            fi.write("gens={}\n\n".format(genes))
            if "copiar" in acciones :
                #Crear la funcion que copia los datos (FASTQs) en la carpeta de analisis
                fi.write("function copiar {\n")
                fi.write("\tcd {}\n".format(pathAnalisi))
                fi.write("\tmkdir {tanda} ; cd {tanda}\n\n".format(tanda = tanda))
                fi.write("\techo -e \"################################\\n\\tCopiant dades\\n################################\\n\"\n")
                for f in fastqs :
                    patx = f.replace(" ", "\\ ") #Parche para leer los espacios en la terminal bash
                    fi.write("\trsync -aP {} .\n".format(patx))

                fi.write("\trsync -aP {} .\n".format(genes))
                fi.write("\tmv ../{} .\n".format(arxiu))
                fi.write("}\n\n")

            fi.write("function analisi {\n")
            fi.write("\tforward=$1\n\treverse=$2\n\treadgroup=$3\n\talias=$4\n")
            fi.write("\techo -e \"################################\\n\\tAnalitzant $alias\\n################################\\n\"\n")
            fi.write("\tmkdir $alias\n")
            fi.write("\tcd $alias\n")
            if "fastqc" in acciones :
                fi.write("\n\t# Control de calidad. FastQC\n")
                fi.write("\tmkdir fastqc # Carpeta donde se guardara el control de calidad\n")
                # Recoger la cadena para invocar fastqc usando la libreria de comandos
                fi.write("\t" + cmd.getFastQC("../$forward", "fastqc") + "\n")
                fi.write("\t" + cmd.getFastQC("../$reverse", "fastqc") + "\n")
                fi.write("\trm fastqc/*zip # Eliminar los archivos comprimidos, ya se han descomprimido al finalizar FastQC\n")
            if "aln" in acciones :
                fi.write("\n\t# Alineamiento. BWA\n")
                # La cadena align tiene cuatro variables: rg es para introducir el read group, fw es para el fastq forward, rv es para el fastq reverse y ref es para el genoma de referencia
                fi.write("\t" + cmd.getAln("$readgroup", "$ref", "../$forward", "../$reverse", "bwa.sam") + "\n")
                fi.write("\t" + cmd.getPcSort("bwa.sam", "bwa.sort.bam") + "\n")
                fi.write("\t" + cmd.getPcIndex("bwa.sort.bam") + "\n")
                fi.write("\tmkdir alignment\n")
                fi.write("\tmv bwa.sam *bam *bai alignment/\n")
            fi.write("\tcd alignment\n")
            # Convertir el bam ordenado en un bed para poder hacer un control de calidad posterior
            if "bamqc" in acciones :
                fi.write("\t" + cmd.getBam2bed("bwa.sort.bam", "bwa.bed") + "\n")
            # Recalibrar las bases
            if "recal" in acciones :
                fi.write("\n\t# Post-Alineamiento. GATK\n")
                # GATKrecal devuelve dos ordenes, separadas por \n. Como queremos que se tabule tambien la segunda linea, se reemplaza el \n por \n + \t
                fi.write("\t" + cmd.getGATKrecal("bwa.sort.bam", "$ref", "$sites", "bwa.recal.bam", "$mani").replace("\n", "\n\t") + "\n")

            # Marcar duplicados
            if "mdups" in acciones :
                fi.write("\t" + cmd.getPcMarkduplicates("bwa.recal.bam", "bwa.nodup.bam") + "\n")
                fi.write("\t" + cmd.getPcIndex("bwa.nodup.bam") + "\n")
            fi.write("\tcd ..\n")
            # Estudios de coverage, on target, off target, porcentaje de bases con X coverage...
            if "coverage" in acciones :
                fi.write("\n\t# Control de calidad del alineamiento y estudio de coverage\n")
                fi.write("\t" + cmd.getCoverageAll("$mani","alignment/bwa.recal.bam", "coverage.txt") + "\n")
                fi.write("\t" + cmd.getCoverageBase("$mani", "alignment/bwa.recal.bam", "coverageBase.txt") + "\n")
                fi.write("\tgrep '^all' coverage.txt > coverageAll.txt\n")
                fi.write("\trm coverage.txt\n")
                fi.write("\tRscript {}/coveragePanells.R\n".format(wd))
                fi.write("\tmkdir coverage\n")
                fi.write("\tmv *png coverageAll.txt coverageBase.txt coverage/\n")
            if "bamqc" in acciones :
                fi.write("\tpython3 {}/bamQC.py\n".format(wd)) # Hay una opcion de lanzar pctDups (calcular porcentaje de duplicados) en caso de exomas
            # Variant calling. La carpeta donde se guardan los datos se llama variantCalling. En caso de queren cambiarse, modificar las dos siguientes lineas
            if "strelkaGerm" in acciones :
                fi.write("\n\t# Variant calling. Strelka2\n")
                if "mdups" in acciones :
                    fi.write("\t" + cmd.getStrelka2("alignment/bwa.nodup.bam", "$ref", "variantCalling", "$gzmani", "variantCalling").replace("\n", "\n\t") + "\n")
                elif "recal" in acciones:
                    fi.write("\t" + cmd.getStrelka2("alignment/bwa.recal.bam", "$ref", "variantCalling", "$gzmani", "variantCalling").replace("\n", "\n\t") + "\n")
                else :
                    fi.write("\t" + cmd.getStrelka2("alignment/bwa.sort.bam", "$ref", "variantCalling", "$gzmani", "variantCalling").replace("\n", "\n\t") + "\n")
                fi.write("\trsync -aP results/variants/variants.vcf.gz .\n")
                fi.write("\tgunzip variants.vcf.gz\n")
                fi.write("\tmv variants.vcf strelka2.vcf")
                # Anotacion de variantes usando ANNOVAR
                if "vanno" in acciones :
                    fi.write("\n\t# Anotacion y filtrado de variantes. ANNOVAR y myvariant.info\n")
                    fi.write("\t" + cmd.getANNOVAR("strelka2.vcf", "raw").replace("\n", "\t\n") + "\n")
                # Re-anotacion y filtrado de variantes usando myvariant.info
                if "filtrar" in acciones :
                    fi.write("\tpython3 {wd}/filtStrelka.py {input} {samplename}".format(wd = wd, input = "raw.hg19_multianno.txt", samplename = "$alias"))
            elif "mutectGerm" in acciones :
                fi.write("\n\t#Variant calling. Mutect2\n")
                if "mdups" in acciones :
                    fi.write("\t" + cmd.getMutect2("alignment/bwa.nodup.bam", "$ref", "$sites", "mutect", "$mani").replace("\n", "\n\t") + "\n")
                elif "recal" in acciones :
                    fi.write("\t" + cmd.getMutect2("alignment/bwa.recal.bam", "$ref", "$sites", "mutect", "$mani").replace("\n", "\n\t") + "\n")
                else :
                    fi.write("\t" + cmd.getMutect2("alignment/bwa.sort.bam", "$ref", "$sites", "mutect", "$mani").replace("\n", "\n\t") + "\n")
                # Anotacion de variantes usando ANNOVAR
                if "vanno" in acciones :
                    fi.write("\n\t# Anotacion y filtrado de variantes. ANNOVAR y myvariant.info\n")
                    fi.write("\t" + cmd.getANNOVAR("mutect.vcf", "raw").replace("\n", "\t\n") + "\n")
                # Re-anotacion y filtrado de variantes usando myvariant.info
                if "filtrar" in acciones :
                    if "mutectGerm" in acciones :
                        fi.write("python3 {wd}/filtMutect.py {input} {samplename}\n".format(wd = wd, input = "raw.hg19_multianno.txt", samplename = "$alias"))


            # Juntar los resultados obtenidos en un Excel
            if "excel" in acciones :
                fi.write("\tpython3 {wd}/data2excel.py {output}\n".format(wd = wd, output = "$alias"))
            fi.write("\tcd {}/{}\n".format(pathAnalisi, tanda))

            fi.write("}\n\n")

            #Crear los comandos para lanzar las funciones
            if "copiar" in acciones :
                fi.write("copiar\n")
            hechos = []
            for f in fastqs :
                params, id = extraerRG(os.path.basename(f))
                if params == "No valido" :
                    print("WARNING: Formato de FASTQ no reconocido para el archivo: {}. Se debera montar la orden para analizar manualmente".format(f))
                else :
                    if id not in hechos :
                        fi.write("analisi {params}\n".format(params = params))
                        hechos.append(id)
        print("INFO: Log guardado como {}".format(arxiu))
        os.chmod(arxiu, 0o754) # Dar permiso de ejecuion para el propietario, lectura y escritura para el grupo, y lectura para el resto de usuarios
