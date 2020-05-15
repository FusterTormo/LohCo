#!/usr/bin/python
# -*- coding: utf-8 -*-

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

def getFastQC(input, output) :
    """Comando para ejecutar FastQC (control de calidad de los FASTQ).

    Los parametros indican
        -o ruta donde se guardaran los archivos de salida.
        -f que el archivo de entrada es un FASTQ
        -extract descomprimir el archivo de salida
        -q omite los mensajes de progreso (log)
        -t 6 el numero de hilos (threads) que usa el programa para ejcutarse en paralelo

    Parameters
    ----------
        input : str
            Ruta al archivo FASTQ que se incrustara en el comando
        output : str
            Carpeta donde se guardara la salida del programa

    Returns
    -------
        str
            Comando para ejecutar FastQC en el archivo que se ha pasado por parametro
    """
    return "/opt/FastQC/fastqc -o {fastqc}/ -f fastq -extract -q -t 6 {fastq}".format(fastq = input, fastqc = output)

def getAln(rg, ref, fw, rv, output) :
    """Comando para ejecutar BWA (alineamiento).

    Los parametros indican
        -M para compatibilidad con Picard tools y GATK
        -t numero de hilos (threads) que usa el programa para ejecutarse
        -R Read Group que se pondra en el sam de salida. Este Read Group es necesario para poder ejecutar GATK (post-alineamiento)

    """
    return "/opt/bwa.kit/bwa mem -M -t 6 -R {rg} {ref} {fw} {rv} > {out}".format(rg = rg, ref = ref, fw = fw, rv = rv, out = output)

def getPcSort(input, output) :
    """Comando para ordenar el bam, usando Picard tools SortBam

    Los parametros indican
        INPUT archivo de entrada que se quiere ordenar
        OUTPUT nombre que tendra el archivo resultante
        SORT_ORDER tipo de ordenamiento que se quiere hacer en el archivo de salida
    """
    return "java -jar /opt/picard-tools-2.21.8/picard.jar SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate".format(input = input, output = output)

def getPcMerge(input, output) :
    """Comando para ordenar una lista de sams

    Los parametros indican
        INPUT cada uno de los sam de entrada que se quieren ordenar
        OUTPUT nombre que tendra el archivo bam resultante
        SORT_ORDER tipo de ordenamiento que se quiere hacer en el archivo de salida
        USE_THREADING flag para lanzar el comando en paralelo
        TMP_DIR directorio donde se guardaran los archivos temporales
    """
    cmd = "java -jar /opt/picard-tools-2.21.8/picard.jar MergeSamFiles "
    for i in input :
        cmd += "INPUT={} ".format(i)
    cmd += "OUTPUT={} SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=/home/ffuster/share/gsole".format(output)
    return cmd

def getPcIndex(bam) :
    """Comando para crear un indice en el bam ordenado
    Los parametros indican
        INPUT archivo de entrada que se quiere indexar
    """
    return "java -jar /opt/picard-tools-2.21.8/picard.jar BuildBamIndex INPUT={bam}".format(bam = bam)

def getBam2bed(bam, bed) :
    """Comando para crear un bed con todas las regiones donde se han alineado reads
    Los parametros indican
        -i archivo (bam) de entrada que se quiere convertir en formato bed
    """
    return "bedtools bamtobed -i {bam} > {bed}".format(bam = bam, bed = bed)

def getGATKrecal(bam, ref, sites, output, regions = None) :
    # Comando para realizar el primer paso de la recalibracion de bases sugerida por GATK
    cmd = "/opt/gatk-4.1.4.1/gatk BaseRecalibrator -I {bam} -R {ref} --known-sites {dbsnp} -O recaldata.table".format(bam = bam, ref = ref, dbsnp = sites)
    if regions != None :
        cmd += " -L {mani}".format(mani = regions)
    cmd += "\n"
    # Comando para realizar el segundo paso de la recalibracion de bases sugerida por GATK
    cmd += "/opt/gatk-4.1.4.1/gatk ApplyBQSR -I {bam} -R {ref} -bqsr-recal-file recaldata.table -O {output}".format(bam = bam, ref = ref, output = output)
    if regions != None :
        cmd += " -L {mani}".format(mani = regions)
    return cmd

def getPcMarkduplicates(input, output) :
    """Comando para marcar duplicados usando Picard tools

    Los parametros indican
        INPUT archivo (bam) en el que se quieren marcar duplicados
        OUTPUT nombre que tendra el bam con los duplicados marcados
        METRICS_FILE archivo donde se muestran los reads que se han marcado como duplicados
    """
    return "java -jar /opt/picard-tools-2.21.8/picard.jar MarkDuplicates INPUT={bamIn} OUTPUT={bamOut} METRICS_FILE=dups_bam.txt".format(bamIn = input, bamOut = output)

def getCoverageAll(regiones, bam, salida) :
    """Comando para calcular el coverage agrupado del bam en las regiones del manifest
    Los parametros indican
        -hist la salida tendra formato histograma
        -a regiones en las que se quiere calcular el coverage
        -b bam con las regiones en las que se calculara el coverage
    """
    return "bedtools coverage -hist -a {mani} -b {bam} > {output}".format(mani = regiones, bam = bam, output = salida)

def getCoverageBase(regiones, bam, salida) :
    """Comando para calcular el coverage de cada una de las bases dentro de la region de interes
    Los parametros indican
        -d la salida mostrara el coverage base por base
        -a regiones en las que se quiere calcular el coverage
        -b bam con las regiones en las que se quiere calcular el coverage

    """
    return "bedtools coverage -d -a {mani} -b {bam} > {output}".format(mani = regiones, bam = bam, output = salida)

def getStrelka2(bam, referencia, output, regiones = None) :
    """ Comando para ejecutar el variant caller Strelka2"""
    cmd = "/opt/strelka-2.9.10/bin/configureStrelkaGermlineWorkflow.py --bam {bam} --referenceFasta {ref} --runDir {variantDir}".format(bam = bam, ref = referencia, variantDir = output)
    if regiones != None :
        cmd += " --exome --callRegions {mani}"
    cmd += "\ncd {variantDir}\n".format(output)
    cmd += "./runWorkflow.py -m local -j 6 --quiet"
    return cmd

def getMutect2(input, referencia, germline, output, regiones = None) :
    cmd = "/opt/gatk-4.1.4.1/gatk Mutect2 -I {input} -R {ref} --germline-resource {gnomad} -O {salida}.vcf".format(input = input, ref = referencia, gnomad = germline, salida = output)
    if regiones != None :
        cmd += " -L {manifest}".format(manifest = regiones)
    cmd += "\n"
    cmd += "/opt/gatk-4.1.4.1/gatk FilterMutectCalls -R {ref} -V {id}.vcf --stats {id}.vcf.stats --filtering-stats {id}.filter.tsv -O {id}.filtered.vcf".format(ref = referencia, id = output)
    return cmd

def getANNOVAR(vcf, prefix) :
    cmd = "/opt/annovar20180416/convert2annovar.pl -format vcf4 -outfile {out}.av -includeinfo {input}\n".format(out = prefix, input = vcf)
    cmd += "/opt/annovar20180416/annotate_variation.pl -geneanno -buildver hg19 -hgvs -separate -out {prefix} {prefix}.av /opt/annovar20180416/humandb/\n".format(prefix = prefix)
    cmd += "/opt/annovar20180416/table_annovar.pl {input}.av /opt/annovar20180416/humandb/ -buildver hg19 -out {input} -remove -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad211_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20190305,cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA -otherinfo".format(input = prefix)
    return cmd
