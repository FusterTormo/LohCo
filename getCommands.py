#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
CONSTANTS:
    Rutas de los programas y parametros de cada uno de los pasos dentro de la pipeline
"""
# IDEA: Podrien ser getters amb els parametres per defecte. Seria mes llegible

picardSort = "java -jar /opt/picard-tools-2.21.8/picard.jar SortSam INPUT=bwa.sam OUTPUT=bwa.sort.bam SORT_ORDER=coordinate" # Comando para ordenar el bam
picardIndex = "java -jar /opt/picard-tools-2.21.8/picard.jar BuildBamIndex INPUT={bam}" # Comando para crear un indice en el bam ordenado
bedtoolsBam2Bed = "bedtools bamtobed -i {bam} > bwa.bed" #Comando para crear un bed con todas las regiones donde se han alineado reads
gatk1 = "/opt/gatk-4.1.4.1/gatk BaseRecalibrator -I {bam} -R {ref} --known-sites {dbsnp} -O recaldata.table" # Comando para realizar el primer paso de la recalibracion de bases sugerida por GATK
gatk2 = "/opt/gatk-4.1.4.1/gatk ApplyBQSR -I {bam} -R {ref} -bqsr-recal-file recaldata.table -O bwa.recal.bam" # Comando para realizar el segundo paso de la recalibracion de bases sugerida por GATK
bedtoolsCoverageAll = "bedtools coverage -hist -a {mani} -b {bam} > {output}" # Comando para calcular el coverage agrupado del bam en las regiones del manifest
bedtoolsCoverageBase = "bedtools coverage -d -a {mani} -b {bam} > {output}" # Comando para calcular el coverage de cada una de las bases dentro de la region de interes
markDup = "java -jar /opt/picard-tools-2.21.8/picard.jar MarkDuplicates INPUT={bam} OUTPUT=bwa.nodup.bam METRICS_FILE=dups_bam.txt" # Comando para marcar duplicados usando Picard tools
vc1 = "/opt/strelka-2.9.10/bin/configureStrelkaGermlineWorkflow.py --bam {bam} --referenceFasta {ref} --exome --runDir {variantDir} --callRegions {mani}" # Comando para ejecutar el variant caller que se va a usar (Strelka2)
vc2 = "./runWorkflow.py -m local -j 6 --quiet"
annovar = "/opt/annovar20180416/convert2annovar.pl -format vcf4 -outfile {out} -includeinfo {vcf}"
annovar2 = "/opt/annovar20180416/annotate_variation.pl -geneanno -buildver hg19 -hgvs -separate -out {output} {input} /opt/annovar20180416/humandb/"
annovar3 = "/opt/annovar20180416/table_annovar.pl {input} /opt/annovar20180416/humandb/ -buildver hg19 -out {output} -remove -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad211_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20190305,cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA -otherinfo"

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

def getFastQC(f) :
    """Comando para ejecutar FastQC (control de calidad de los FASTQ).

    Los parametros indican
        -o ruta donde se guardaran los archivos de salida.
        -f que el archivo de entrada es un FASTQ
        -extract descomprimir el archivo de salida
        -q omite los mensajes de progreso (log)
        -t 6 el numero de hilos (threads) que usa el programa para ejcutarse en paralelo

    Parameters
    ----------
        f : str
            Ruta al archivo FASTQ que se incrustara en el comando

    Returns
    -------
        str
            Comando para ejecutar FastQC en el archivo que se ha pasado por parametro
    """
    return "/opt/FastQC/fastqc -o fastqc/ -f fastq -extract -q -t 6 {fastq}".format(f)

def getAln(rg, ref, fw, rv) :
    """Comando para ejecutar BWA (alineamiento).

    Los parametros indican
        -M para compatibilidad con Picard tools y GATK
        -t numero de hilos (threads) que usa el programa para ejecutarse
        -R Read Group que se pondra en el sam de salida. Este Read Group es necesario para poder ejecutar GATK (post-alineamiento)

    """
    return "/opt/bwa.kit/bwa mem -M -t 6 -R {rg} {ref} {fw} {rv} > bwa.sam".format(rg = rg, ref = ref, fw = fw, rv = rv)
