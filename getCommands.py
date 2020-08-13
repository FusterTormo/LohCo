#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
CONSTANTS:
    Rutas de los programas de cada uno de los pasos dentro de la pipeline
"""
fastqc = "/opt/FastQC/fastqc"
bwa = "/opt/bwa.kit/bwa"
picard = "/opt/picard-tools-2.21.8/picard.jar"
bedtools = "bedtools"
gatk = "/opt/gatk-4.1.8.0/gatk"
strelka2 = "/opt/strelka-2.9.10"
samtools = "samtools"
varscan = "/opt/varscan-2.4.0/VarScan.v2.4.0.jar"
annovar = "/opt/annovar20200607"

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
    return "{exec} -o {fastqc}/ -f fastq -extract -q -t 6 {fastq}".format(fastq = input, fastqc = output, exec = fastqc)

def getAln(rg, ref, fw, rv, output) :
    """Comando para ejecutar BWA (alineamiento).

    Los parametros indican
        -M para compatibilidad con Picard tools y GATK
        -t numero de hilos (threads) que usa el programa para ejecutarse
        -R Read Group que se pondra en el sam de salida. Este Read Group es necesario para poder ejecutar GATK (post-alineamiento)

    Parameters
    ----------
        rg : str
            Read group. Formato minimo para ser aceptado por GATK: "@RG\tID:___\tSM:____\tPL:____"
        ref : str
            Ruta donde esta el genoma de referencia
        fw : str
            Ruta donde esta el fastq forward (o read 1)
        rv : str
            Ruta donde esta el fastq reverse (o read 2)
        output : str
            Ruta y nombre donde se guardara el sam resultante del alineamiento

    Returns
    -------
        str
            Comando para ejecutar BWA en los fastq que se han pasado por parametro
    """
    return "{exec} mem -M -t 6 -R {rg} {ref} {fw} {rv} > {out}".format(rg = rg, ref = ref, fw = fw, rv = rv, out = output, exec = bwa)

def getPcSort(input, output) :
    """Comando para ordenar el bam, usando Picard tools SortBam

    Los parametros indican
        INPUT archivo de entrada que se quiere ordenar
        OUTPUT nombre que tendra el archivo resultante
        SORT_ORDER tipo de ordenamiento que se quiere hacer en el archivo de salida

    Parameters
    ----------
        input : str
            Ruta del sam que se quiere ordenar
        output : str
            Nombre y ruta que tendra el bam ordenado resultante

    Returns
    -------
        str
            Comando para ejecutar Picard SortBam en el archivo bam que se ha pasado por parametro
    """
    return "java -jar {picard} SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate".format(input = input, output = output, picard = picard)

def getPcMerge(input, output) :
    """Comando para ordenar y agrupar una lista de sams

    Los parametros indican
        INPUT cada uno de los sam de entrada que se quieren ordenar
        OUTPUT nombre que tendra el archivo bam resultante
        SORT_ORDER tipo de ordenamiento que se quiere hacer en el archivo de salida
        USE_THREADING flag para lanzar el comando en paralelo
        TMP_DIR directorio donde se guardaran los archivos temporales

    Parameters
    ----------
        input : list
            Lista de rutas donde estan cada uno de los sams que se quiere agrupar/ordenar
        output : str
            Nombre y ruta que tendra el bam ordenado resultante

    Returns
    -------
        str
            Comando para ejecutar Picard MergeSamFiles la lista de sams que se ha pasado por parametro
    """
    cmd = "java -jar {exec} MergeSamFiles ".format(exec = picard)
    for i in input :
        cmd += "INPUT={} ".format(i)
    cmd += "OUTPUT={} SORT_ORDER=coordinate USE_THREADING=true TMP_DIR=/home/ffuster/share/gsole".format(output)
    return cmd

def getPcIndex(bam) :
    """Comando para crear un indice en un bam ordenado

    Los parametros indican
        INPUT archivo de entrada que se quiere indexar

    Parameters
    ----------
        bam : str
            Ruta donde se encuentra el bam al que se quiere crear el indice. El archivo resultante se guardara en el mismo sitio, pero con extension .bai en lugar de .bam

    Returns
    -------
        str
            Comando para ejecutar Picard BuildBamIndex en el archivo que se ha pasado por parametro
    """
    return "java -jar {picard} BuildBamIndex INPUT={bam}".format(bam = bam, picard = picard)

def getBam2bed(bam, bed) :
    """Comando para crear un bed con todas las regiones donde se han alineado reads

    Los parametros indican
        -i archivo (bam) de entrada que se quiere convertir en formato bed

    Parameters
    ----------
        bam : str
            Ruta del archivo bam en el que se quiere crear el bed con sus regiones
        bed : str
            Nombre y ruta del bed resultante

    Returns
    -------
        str
            Comando para ejecutar bedtools bamtobed en el archivo que se ha pasado por parametro
    """
    return "{exec} bamtobed -i {bam} > {bed}".format(bam = bam, bed = bed, exec = bedtools)

def getGATKrecal(bam, ref, sites, output, regions = None) :
    """Comando para realizar el primer paso de la recalibracion de bases sugerida por GATK. Este comando se divide en dos: BaseRecalibrator i ApplyBQSR

    BaseRecalibrator genera el recalibrado basandose en distintas covariaciones. Los parametros indican
        -I bam que se quiere recalibrar
        -R genoma de referencia en el que se ha alineado dicho bam
        --known-sites vcf con SNPs recogidas en bases de datos poblacionales. Se pueden descargar desde https://gnomad.broadinstitute.org/downloads
        -O tabla donde se guardaran los datos de recalibrado para su posterior uso en ApplyBQSR
    ApplyBQSR aplica el recalibrado de bases en el model que se ha entrenado al ejecutar BaseRecalibrator. Los parametros indican
        -I bam que se quiere recalibrar
        -R genoma de referencia en el que se ha alineado el bam
        --bqsr-recal-file tabla de recalibrado que ha generado BaseRecalibrator
        -O nombre del bam de salida que generara el programa. Ademas del bam, el programa crea un indice

    Parameters
    ----------
        bam : str
            Ruta donde esta el bam que se quiere recalibrar
        ref : str
            Ruta donde esta el genoma de referencia en el que se ha alineado dicho bam
        sites : str
            Ruta donde esta el vcf con los snps que se usara para en BaseRecalibrator
        output : str
            Nombre y ruta que tendra el bam recalibrado que generen esos comandos
        regions : str, optional
            Ruta donde esta el archivo bed con las regiones en las que se quiere focalizar el recalibrado. En caso de paneles se puede usar el manifest file

    Returns
    -------
        str
            Comando para ejecutar el recalibrado de bases sugerido por GATK en el archivo bam que se ha pasado por parametro
    """
    cmd = "{gatk} BaseRecalibrator -I {bam} -R {ref} --known-sites {dbsnp} -O recaldata.table".format(bam = bam, ref = ref, dbsnp = sites, gatk = gatk)
    if regions != None :
        cmd += " -L {mani}".format(mani = regions)
    cmd += "\n"
    # Comando para realizar el segundo paso de la recalibracion de bases sugerida por GATK
    cmd += "{gatk} ApplyBQSR -I {bam} -R {ref} -bqsr-recal-file recaldata.table -O {output}".format(bam = bam, ref = ref, output = output, gatk = gatk)
    if regions != None :
        cmd += " -L {mani}".format(mani = regions)
    return cmd

def getPcMarkduplicates(input, output) :
    """Comando para marcar duplicados usando Picard tools

    Los parametros indican
        INPUT archivo (bam) en el que se quieren marcar duplicados
        OUTPUT nombre que tendra el bam con los duplicados marcados
        METRICS_FILE archivo donde se muestran los reads que se han marcado como duplicados

    Parameters
    ----------
        input : str
            Ruta donde esta el bam en el que se quieren marcar duplicados
        output : str
            Nombre y ruta del archivo con los duplicados marcados que se generara al ejecutar este comando

    Returns
    -------
        str
            Comando para ejecutar Picard MarkDuplicates en el archivo que se ha pasado por parametro
    """
    return "java -jar {picard} MarkDuplicates INPUT={bamIn} OUTPUT={bamOut} METRICS_FILE=dups_bam.txt".format(bamIn = input, bamOut = output, picard = picard)

def getCoverageAll(regiones, bam, salida) :
    """Comando para calcular el coverage agrupado del bam en las regiones del manifest

    Los parametros indican
        -hist la salida tendra formato histograma
        -a regiones en las que se quiere calcular el coverage
        -b bam con las regiones en las que se calculara el coverage

    Parameters
    ----------
        regiones : str
            Ruta donde esta el archivo bed con las coordenadas donde queremos calcular el coverage
        bam : str
            Ruta donde esta el bam en el que queremos calcular el coverage
        salida : str
            Nombre y ruta del archivo de coverage que se generara al ejecutar este comando

    Returns
    -------
        str
            Comando para ejecutar bedtools coverage en el archivo bam que se ha pasado por parametro
    """
    return "{exec} coverage -hist -a {mani} -b {bam} > {output}".format(mani = regiones, bam = bam, output = salida, exec = bedtools)

def getCoverageBase(regiones, bam, salida) :
    """Comando para calcular el coverage de cada una de las bases dentro de la region de interes

    Los parametros indican
        -d la salida mostrara el coverage base por base
        -a regiones en las que se quiere calcular el coverage
        -b bam con las regiones en las que se quiere calcular el coverage

    Parameters
    ----------
        regiones : str
            Ruta donde esta el archivo bed con las coordenadas donde queremos calcular el coverage
        bam : str
            Ruta donde esta el bam en el que queremos calcular el coverage
        salida : str
            Nombre y ruta del archivo de coverage que se generara al ejecutar este comando

    Returns
    -------
        str
            Comando para ejecutar bedtools coverage en el archivo bam que se ha pasado por parametro
    """
    return "{exec} coverage -d -a {mani} -b {bam} > {output}".format(mani = regiones, bam = bam, output = salida, exec = bedtools)

def getStrelka2(bam, referencia, output, regiones = None) :
    """Comandos para ejecutar el variant caller Strelka2

    El variant calling de Strelka2 se compone de dos pasos. Configuracion y ejecucion.
    En la configuracion se crean los scripts necesarios para poder ejecutar Strelka2 en modo de analisis germinal, en este caso. Los parametros indican
        --bam archivo bam en el que se va a ejecutar el variant calling
        --referenceFasta genoma de referencia en el que se ha alineado el bam
        --runDir carpeta donde se van a guardar todos los scripts necesarios para el analisis. En esta carpeta se guardaran los archivos resultantes del variant calling. En caso de no existir, se crea la carpeta
        --exome indica que el bam no es un genoma (whole genome sequencing), sino un exoma (whole exome sequencing) o un panel (targeted deep sequencing)
        --callRegions regiones en las que se esta interesado en hacer el variant calling. En caso de analisis de paneles, en este parametro se adjuntara el archivo manifest
    En la ejecucion se ejecuta el variant calling propiamente dicho. Los parametro indican
        -m tipo de computadora donde se ejecutara el variant calling. Opciones disponibles: local, sge. La opcion sge suele dar errores de ejecucion incluso cuando se ejecuta en un cluster
        -j numero de hilos (threads, procesos) que usara Strelka2 para hacer el variant calling
        --quiet eliminar los mensajes de log que no sean errores de la terminal

    Parameters
    ----------
        bam : str
            Ruta del archivo bam en el que se quiere hacer el variant calling
        referencia : str
            Genoma de referencia en el que se ha hecho el alineamiento
        output : str
            Nombre de la carpeta donde se guardaran los archivos de ejecucion y resultantes. Si esta carpeta no existe, Sterlka2 la creara
        regiones : str, optional
            Ruta del archivo bed con las coordenadas genomicas donde se quiere hacer el variant calling. Este parametro no se debe pasar a no ser que se esten analizando paneles

    Returns
    -------
        str
            Comandos para ejecutar Strelka2 en el archivo bam que se ha pasado por parametro
    """
    cmd = "{strelka}/bin/configureStrelkaGermlineWorkflow.py --bam {bam} --referenceFasta {ref} --runDir {variantDir}".format(bam = bam, ref = referencia, variantDir = output, strelka = strelka2)
    if regiones != None :
        cmd += " --exome --callRegions {mani}".format(mani = regiones)
    elif regiones == "exoma" :
        cmd += " --exome"
    cmd += "\ncd {variantDir}\n".format(output)
    cmd += "./runWorkflow.py -m local -j 6 --quiet"
    return cmd

def getStrelka2soma(tumor, control, referencia, output, regiones = None) :
    """Comandos para ejecutar el variant caller Strelka2 para analisis somatico

    El variant calling de Strelka2 se compone de 2 pasos. Configuracion y ejecucion.
    En la configuracion se crean los scripts necesarios para poder ejecutar Strelka2 en modo de analisis somatico, en este caso. Los parametros indican
        --tumorBam archivo bam de la muestra tumoral
        --normalBam archivo bam de la muestra control/normal
        --referenceFasta genoma de referencia en el que se han alineado ambos bam
        --runDir carpeta donde se guardaran los scripts necesarios para el analisis. Tambien se guardaran los archivos resultantes del variant calling. Se creara la carpeta si no existe previamente
        --exome indica al programa que el bam no es un genoma (whole genome sequencing), sino un exoma (whole exome sequencing) o un panel (target deep sequencing)
        --callRegions regiones en las que se esta interesado hacer el variant calling. En caso de analisis de paneles, en este parametro de debe adjuntar el archivo manifest
    En la ejecucion se lanza el variant calling propiamente dicho. Los parametros indican
        -m tipo de computadora donde se ejecutara el variant calling. Opciones disponibles: local, sge. La opcion sge suele dar errores de ejecucion incluso cuando se ejecuta en un cluster
        -j numero de hilos (procesos, threads) que usara Strelka2 para hacer el variant calling
        --quiet eliminar los mensajes de log de la consola, a no ser que sean mensajes de error

    Parameters
    ----------
        tumor : str
            Ruta del bam creado para la muestra somatica
        control : str
            Ruta del bam creado para la muestra control (normal)
        referencia : str
            Ruta donde esta el genoma de referencia sobre el que se han alineado las muestras tumoral y normal
        output : str
            Nombre y ruta de la carpeta donde se guardaran los archivos resultantes y scripts generados por Strelka2. Si esta carpeta no existe, se creara automaticamente
        regiones : str, optional
            Ruta del archivo bed con las coordenadas genomicas donde se quiere hacer el variant calling. Este parametro no se debe pasar a no ser que se esten analizando paneles

    Returns
    -------
        str
            Comandos para ejecutar Strelka2 en el archivo bam que se ha pasado por parametro
    """
    cmd = "{strelka}/bin/configureStrelkaSomaticWorkflow.py --tumorBam {tm} --normalBam {cn} --referenceFasta {ref} --runDir {varDir}".format(tm = tumor, cn = control, ref = referencia, varDir = output, strelka = strelka2)
    if regiones != None :
        cmd += " --exome --callRegions{mani}".format(mani = regiones)
    elif regiones == "exoma" :
        cmd += " --exome"
    cmd += "\ncd {variantDir}\n".format(variantDir = output)
    cmd += "./runWorkflow.py -m local -j 6 --quiet"
    return cmd

def getMutect2(input, referencia, germline, output, regiones = None) :
    """Comandos para ejecutar el variant caller Mutect2 en analisis somatico, pero sin muestra control

    El variant calling usando Mutect2 se compone de dos pasos. Variant calling y filtrado a posteriori.
    El paso de variant calling hace la llamada de variantes somaticas, en este caso, sin tener un control asociado. Este analisis es parecido al analisis germinal usado por Strelka2 y VarScan2.
    Los parametros indican
        -I archivo bam en el que se quiere hacer el variant calling
        -R genoma de referencia que se ha usado para generar el bam
        --germline-resource vcf con SNPs recogidas de bases de datos poblacionales
        -O nombre del archivo de salida
    GATK ha creado, ademas, un programa para filtrar las variantes una vez se han llamado. Los parametros indican
        -R genoma de referencia que se ha usado para generar el bam
        -V nombre del archivo vcf que ha creado Mutect2
        --stats nombre del archivo de estadisticas que ha generado Mutect2 a partir del vcf
        --filtering-stats nombre del archivo de estadisticas de filtros que ha generado Mutect2 a partir del vcf
        -O nombre del archivo vcf filtrado que se generara al final del programa

    Parameters
    ----------
        input : str
            Ruta del archivo bam en el que se quiere hacer el variant calling
        referencia : str
            Ruta del genoma de referencia usado para crear el bam que se pasa como parametro
        germline : str
            Ruta del vcf de SNPs recogidas de las bases de datos poblacionales que se usaran para llamar variantes
        output : str
            Nombre y ruta donde se guardaran los archivos resultantes de la ejecucion del variant calling
        regiones : str, optional
            Regiones de interes para el variant calling. Solo se deberia usar en caso de que la muestra a analizar sea un panel de secuenciacion

    Returns
    -------
        str
            Comandos para ejecutar Mutect2 en el archivo bam que se ha pasado por parametro
    """
    cmd = "{gatk} Mutect2 -I {input} -R {ref} --germline-resource {gnomad} -O {salida}.vcf".format(input = input, ref = referencia, gnomad = germline, salida = output, gatk = gatk)
    if regiones != None :
        cmd += " -L {manifest}".format(manifest = regiones)
    cmd += "\n"
    cmd += "{gatk} FilterMutectCalls -R {ref} -V {id}.vcf --stats {id}.vcf.stats --filtering-stats {id}.filter.tsv -O {id}.filtered.vcf\n".format(ref = referencia, id = output, gatk = gatk)
    cmd += "mkdir variantCalling\n"
    cmd += "cd variantCalling\n"
    cmd += "mv ../mutect* ."
    return cmd

def getMutect2soma(tumor, control, idtum, idcon, referencia, germline, outputdir, regiones = None) :
    """Comandos para ejecutar el variant caller Mutect2 en analisis somatico tradicional

    El paso de variant calling hace la llamada de variantes somaticas, en este caso, sin tener un control asociado. Este analisis es parecido al analisis germinal usado por Strelka2 y VarScan2.
    Los parametros indican
        --java-options indica opciones adicionales al interprete de java. En este caso se reservan 32Gb de memoria RAM para la ejecucion del programa
        -I archivos bam en los que se va a ejecutar el variant calling. Mutect2 considera tanto la muestra control como la tumoral como bam de entrada
        -tumor nombre de la muestra, especificado en el read group del alineamiento, que se debe considerar como muestra tumoral. Para obtener este dato, se puede ejecutar samtools view -H {tumor}
            y buscar dentro de la cabecera la etiqueta @RG (continene los datos del read group). Luego buscar el contenido de la etiqueta SM
        -normal nombre de la muestra, especificado en el read group del alineamiento, que se debe considerar como muestra normal
        -R genoma de referencia que se ha usado para generar el bam
        --germline-resource vcf con SNPs recogidas de bases de datos poblacionales
        -O nombre del archivo de salida
    GATK ha creado, ademas, un programa para filtrar las variantes una vez se han llamado. Los parametros indican
        -R genoma de referencia que se ha usado para generar el bam
        -V nombre del archivo vcf que ha creado Mutect2
        --stats nombre del archivo de estadisticas que ha generado Mutect2 a partir del vcf
        --filtering-stats nombre del archivo de estadisticas de filtros que ha generado Mutect2 a partir del vcf
        -O nombre del archivo vcf filtrado que se generara al final del programa

    Parameters
    ----------
        tumor : str
            Ruta del archivo bam en el que se quiere hacer el variant calling. Este bam debe ser el de la muestra tumoral
        control : str
            Ruta del archivo bam en el que se quiere hacer el variant calling. Este bam debe ser el de la muestra control
        idtum : str
            Valor de la etiqueta SM de la muestra tumoral
        idcon : str
            Valor de la etiqueta SM de la muestra normal
        referencia : str
            Ruta del genoma de referencia usado para crear el bam que se pasa como parametro
        germline : str
            Ruta del vcf de SNPs recogidas de las bases de datos poblacionales que se usaran para llamar variantes
        outputdir : str
            Nombre y ruta de la carpeta donde se guardaran los archivos resultantes de la ejecucion del variant calling. La creacion de esta carpeta viene adjuntada en los comandos para lanzar Mutect2
        regiones : str, optional
            Regiones de interes para el variant calling. Solo se deberia usar en caso de que la muestra a analizar sea un panel de secuenciacion

    Returns
    -------
        str
            Comandos para ejecutar Mutect2 en los archivos bam que se han pasado por parametro
    """
    cmd = "{gatk} --java-options \"-Xmx32G\" Mutect2 -I {tm} -I {cn} -tumor {idtm} -normal {idcn} -R {ref} --germline-resource {gnomad} -O {idtm}.vcf".format(tm = tumor, cn = control, idtm = idtum, idcn = idcon, ref = referencia, gnomad = germline, gatk = gatk)
    if regiones != None :
        cmd += " -L {manifest}".format(manifest = regiones)
    cmd += "\n"
    cmd += "{gatk} FilterMutectCalls -R {ref} -V {id}.vcf --stats {id}.vcf.stats --filtering-stats {id}.filter.tsv -O {id}.filtered.vcf\n".format(ref = referencia, id = idtum, gatk = gatk)
    cmd += "mkdir {dir}\n".format(dir = outputdir)
    cmd += "mv {idtm}* {dir}".format(idtm = idtum, dir = outputdir)
    return cmd

def getVarScan2Soma(tumor, control, referencia, output = "VarScan2") :
    """Comandos para ejecutar el variant caller VarScan2 para analisis somatico

    Para ejecutar VarScan2, se tiene que crear un archivo pileup previamente. Este archivo se crea ejecutando samtools mpileup. Los parametros indican
        -f genoma de referencia que se ha usando para crear el bam en el que se va a hacer el variant calling
        -q calidad minima de alineamiento que debe tener la variante para que sea considerada
    Una vez creados los pileup, se ejecuta VarScan2. Notese que VarScan2 se compila al vuelo, no tiene un ejecutable propiamente dicho. Los parametros indican
        --output-vcf indica que el formato de la salida debe ser vcf estandar
        --strand-filter indica a VarScan2 que debe usar el filtro de strand-bias para la llamada de variantes (variant calling)

    Parameters
    ----------
        tumor : str
            Ruta donde esta el bam de la muestra tumoral que se va a usar para el variant calling
        control : str
            Ruta donde esta el bam de la muestra control que se va a usar para el variant calling
        referencia : str
            Ruta donde se encuentra el genoma de referencia que se ha usado para alinear los bams
        output : str, optional
            Nombre de la carpeta donde se guardaran los archivos resultantes del variant calling. Por defecto "VarScan2"

    Returns
    -------
        str
            Comandos para ejecutar VarScan2 en los archivos bam que se han pasado por parametro
    """
    cmd = "{samtools} mpileup -f {ref} -q 1 {bam} > control.pileup\n".format(ref = referencia, bam = control, samtools = samtools)
    cmd += "{samtools} mpileup -f {ref} -q 1 {bam} > tumor.pileup\n".format(ref = referencia, bam = tumor, samtools = samtools)
    cmd += "java -jar {varscan} somatic control.pileup tumor.pileup varscan --output-vcf 1 --strand-filter 1\n".format(varscan = varscan)
    cmd += "mkdir {dir}\n".format(dir = output)
    cmd += "mv tumor.pileup control.pileup varscan* {dir}".format(dir = output)
    return cmd

def getANNOVAR(vcf, prefix, hgref = "hg19", format = "vcf4") :
    """Comandos para anotar las variantes usando ANNOVAR

    ANNOVAR usa dos pasos para anotar las variantes. Convertir el vcf a formato ANNOVAR y anotar dicho archivo
    La conversion transforma el archivo vcf en un formato que sea legible para ANNOVAR posteriormente. Los parametros indican
        -format tipo de vcf que se ha creado y que se quiere convertir a ANNOVAR. Los formatos mas tipicos son vcf4 (muestras germinales) y vcf4old (en caso de vcfs de muestras somaticas)
        -outfile nombre que se quiere dar al archivo convertido a formato ANNOVAR
        -includeinfo incluir en el archivo con formato ANNOVAR la infomacion completa del vcf original
    Una vez se ha convertido el vcf en formato ANNOVAR se ejecuta la anotacion de variantes. Este paso recoge informacion de distintas bases de datos genomicas (clinicas, poblacionales, geneticas...).
        Los parametros indican
        -buildver tipo de genoma de referencia en el que quieren hacer las anotaciones. Se usa para escoger las librerias de anotaciones
        -out prefijo de los archivos de anotacion que se crearan al anotar las variantes
        -remove indica al programa que tiene que eliminar los archivos temporales una vez termine la anotacion
        -protocol nombre de las bases de datos de las que se quiere extraer la informacion de las variantes
        -operation tipo de anotacion que se hara: g indica genomico (posicion que ocupa la coordenada dentro del genoma), f indica filtro, aunque en este caso no se filtre informacion
        -nastring en caso de que no haya informacion de la variante en un campo en concreto de la anotacion, cadena que se escribira
        -otherinfo incluir la informacion del vcf en el archivo de texto final

    Parameters
    ----------
        vcf : str
            Ruta donde esta el archivo vcf que se quiere anotar
        prefix : str
            Nombre que tendran los archivos que genere ANNOVAR
        hgref : str, optional
            Genoma de referencia usado para hacer el analisis del panel. Se para saber de que base de datos recoger la informacion
        format : str, optional
            Tipo de vcf que se ha generado. Este parametro solo hace falta modificarlo en caso de analisis somatico. En dicho caso, se cambiara vcf4 por vcf4old
    """
    if hgref == "hg19" :
        cmd = "{annovar}/convert2annovar.pl -format {format} -outfile {out}.av -includeinfo {input}\n".format(out = prefix, input = vcf, format = format, annovar = annovar)
        cmd += "{annovar}/table_annovar.pl {input}.av /home/ffuster/share/biodata/Indexes/ANNOVAR/humandb/ -buildver hg19 -out {input} -remove -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad211_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20190305,cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA -otherinfo".format(input = prefix, annovar = annovar)
    elif hgref == "hg38" :
        cmd = "{annovar}/convert2annovar.pl -format {format} -outfile {out}.av -includeinfo {input}\n".format(out = prefix, input = vcf, format = format, annovar = annovar)
        cmd += "{annovar}/table_annovar.pl {input}.av /home/ffuster/share/biodata/Indexes/ANNOVAR/humandb/ -buildver hg38 -out {input} -remove -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,gnomad211_exome,gnomad30_genome,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,clinvar_20200316,cosmic70,dbnsfp35a --operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring NA -otherinfo".format(input = prefix, annovar = annovar)
    return cmd
