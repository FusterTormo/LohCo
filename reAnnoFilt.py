#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Re-anotar y filtrar el archivo de salida de table_annovar
"""
import sys

def addFORMAT(dic, somatic) :
    """
    Separar los datos de las columnas INFO y FORMAT del vcf en columnas. Se devuelven los datos separados en un diccionario. Se comprueba si alguna clave esta ya presente en el diccionario
    """
    pass # Recogida desde strelka-allinfo2excel

def convertirData(path) :
    """Abrir un archivo con formato table_annovar y guardar todo el contenido en un diccionario"""
    cabecera = True
    # Columnas dentro del table_annovar
    claves = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'avsnp150', '1000g2015aug_all', '1000g2015aug_afr',
    '1000g2015aug_amr', '1000g2015aug_eas', '1000g2015aug_eur', '1000g2015aug_sas', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS',
    'gnomad_exome_AF', 'gnomad_exome_AF_popmax', 'gnomad_exome_AF_male', 'gnomad_exome_AF_female', 'gnomad_exome_AF_raw', 'gnomad_exome_AF_afr', 'gnomad_exome_AF_sas', 'gnomad_exome_AF_amr',
    'gnomad_exome_AF_eas', 'gnomad_exome_AF_nfe', 'gnomad_exome_AF_fin', 'gnomad_exome_AF_asj', 'gnomad_exome_AF_oth', 'gnomad_exome_non_topmed_AF_popmax', 'gnomad_exome_non_neuro_AF_popmax',
    'gnomad_exome_non_cancer_AF_popmax', 'gnomad_exome_controls_AF_popmax',
    'gnomad_genome_AF', 'gnomad_genome_AF_popmax', 'gnomad_genome_AF_male', 'gnomad_genome_AF_female', 'gnomad_genome_AF_raw', 'gnomad_genome_AF_afr', 'gnomad_genome_AF_sas', 'gnomad_genome_AF_amr',
    'gnomad_genome_AF_eas', 'gnomad_genome_AF_nfe', 'gnomad_genome_AF_fin', 'gnomad_genome_AF_asj', 'gnomad_genome_AF_oth', 'gnomad_genome_non_topmed_AF_popmax', 'gnomad_genome_non_neuro_AF_popmax',
    'gnomad_genome_non_cancer_AF_popmax', 'controls_AF_popmax', 'esp6500siv2_all', 'esp6500siv2_ea', 'esp6500siv2_aa', 'CLNALLELEID', 'CLNDN', 'CLNDISDB',
    'CLNREVSTAT', 'CLNSIG', 'cosmic70', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score',
    'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred',
    'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore',
    'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore',
    'M-CAP_pred', 'REVEL_score', 'REVEL_rankscore', 'MutPred_score', 'MutPred_rankscore', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score',
    'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score',
    'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian',
    'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds',
    'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue', 'CHROM', 'POS', 'ID', 'REFERENCE', 'ALTERATED', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    with open(path, "r") as fi :
        for l in fi :
            if cabecera :
                print(l)
                break
            else :
                pass


def addWebInfo(dc) :
    """Agregar la informacion de myvariant.info a aquellas variantes que sean exonicas"""
    # TODO: Anadir estadisticas de los tipos de variantes recogidos en el panel a la pestaña QC
    pass # Recodigda de strelka2excel


def calcularStrandBias(l, somatic) :
    coef = "NP"
    if somatic == "germinal" :
        auxF = l["ADF"].split(",")
        auxR = l["ADR"].split(",")
        coef = int(auxF[0]) + int(auxF[1]) - int(auxR[0]) - int(auxR[1])
    elif somatic == "somatico" and "SNVSB_INFO" in l.keys() :
        coef = l["SNVSB_INFO"]
    # Las indels somaticas no tienen ninguna forma de calcular el strand bias
    return str(coef)

def getGenomeType(l, somatic) :
    """
    Devolver el tipo de variante que ha detectado el vcf. En caso de que la muestra a anotar sea somatica, se devolvera siempre 0/1 (heterozigoto). Si la muestra es germinal, se devuelve
    el valor de la columna GT del vcf.
    """
    ret = "NP"
    if somatic == "somatico" :
        ret = "0/1"
    else :
        ret = l["GT"]
    return ret

def getGenomeQuality(l, somatic) :
    """
    Devolver la calidad de la variante que ha detectado el variant caller. Strelka2 devuelve la calidad genomica en caso de variant calling somatico y la calidad de la variante en caso de variantes
    somaticas. El nombre de la columna de la variante somatica depende de si la variante es una SNV o una indel.
    """
    ret = "NP"
    if somatic == "germinal" :
        ret = l["GQX"]
    elif somatic == "somatico" and "QSS" in l.keys() : #Calidad de la variante en SNVs somaticas
        ret = l["QSS_INFO"]
    elif somatic == "somatico" and "QSI" in l.keys() : #Calidad de la variante en indels somaticas
        ret = l["QSI_INFO"]
    return ret

def getRegionDepth(l, somatic) :
    """
    Devolver la profundidad (coverage) en la region donde esta la variante. Para variantes germinales, Strelka2 devuelve el dato en las columnas DP, DPF y AD. Como el dato de AD separa las
    coberturas entre referencia y alterado, uso esta columna para calcular la cobertura total. Para variantes somaticas, uso la suma de las coberturas de cada una de las bases, en caso de SNV y
    la suma de las coberturas de referencia y alteradas en caso de indels
    """
    dp = "NP"
    if somatic == "germinal" :
        aux = l["AD"].split(",")
        dp = int(aux[0]) + int(aux[1])
    elif somatic == "somatico" and "AU" in l.keys() :
        dp = int(l["AU_TUMOR"].split(",")[0]) + int(l["CU_TUMOR"].split(",")[0]) + int(l["GU_TUMOR"].split(",")[0]) + int(l["TU_TUMOR"].split(",")[0])
    elif somatic == "somatico" and "TAR" in l.keys() :
        dp = int(l["TAR_TUMOR"].split(",")[0]) + int(l["TIR_TUMOR"].split(",")[0])

    return str(dp)

def getReferenceDepth(l, somatic) :
    """
    Devolver la profundidad (coverage) de reads con la base de referencia. Uso los mismos datos que en getRegionDepth para ahorrarme inconsistencias en el calculo de datos
    """
    dp = "NP"
    if somatic == "germinal" :
        dp = l["AD"].split(",")[0]
    elif somatic == "somatico" and "AU_TUMOR" in l.keys() : #Comprobar si la variante es una SNV somatica
        if l["Ref"] == "A" :
            dp = l["AU_TUMOR"].split(",")[0]
        elif l["Ref"] == "C" :
            dp = l["CU_TUMOR"].split(",")[0]
        elif l["Ref"] == "G" :
            dp = l["GU_TUMOR"].split(",")[0]
        elif l["Ref"] == "T" :
            dp = l["TU_TUMOR"].split(",")[0]
    elif somatic == "somatico" and "TAR_TUMOR" in l.keys() : #Comprobar si la variante es una indel somatica
        dp = l["TAR_TUMOR"].split(",")[0]

    return dp

def getAlteratedDepth(l, somatic) :
    """
    Devolver la profundidad (coverage) de reads con la base alterada. Uso los mismos datos que en getRegionDepth y en getReferenceDepth para ahorrarme inconsistencias en el calculo de datos
    """
    dp = "NP"
    if somatic == "germinal" :
        dp = l["AD"].split(",")[1]
    elif somatic == "somatico" and "AU_TUMOR" in l.keys() : #Comprobar si la variante es una SNV somatica
        if l["Alt"] == "A" :
            dp = l["AU_TUMOR"].split(",")[0]
        elif l["Alt"] == "C" :
            dp = l["CU_TUMOR"].split(",")[0]
        elif l["Alt"] == "G" :
            dp = l["GU_TUMOR"].split(",")[0]
        elif l["Alt"] == "T" :
            dp = l["TU_TUMOR"].split(",")[0]
    elif somatic == "somatico" and "TIR_TUMOR" in l.keys() : #Comprobar si la variante es una indel somatica
        dp = l["TIR_TUMOR"].split(",")[0]

    return dp

def getVAF(l, somatic) :
    """
    Calcula la Variant Allele Frequency VAF usando los datos de coverages de la forma en que se recogen en las funciones anteriores (getRegionDepth, getReferenceDepth, getAlteratedDepth)
    """
    dp = "NP"
    try :
        if somatic == "germinal" :
            aux = l["AD"].split(",")
            dp = float(aux[1]) / (float(aux[0]) + float(aux[1])) * 100
        elif somatic == "somatico" and "AU" in l.keys() :
            ref = 0
            alt = 0
            if l["Ref"] == "A" :
                ref = l["AU_TUMOR"].split(",")[0]
            elif l["Ref"] == "C" :
                ref = l["CU_TUMOR"].split(",")[0]
            elif l["Ref"] == "G" :
                ref = l["GU_TUMOR"].split(",")[0]
            elif l["Ref"] == "T" :
                ref = l["TU_TUMOR"].split(",")[0]

            if l["Alt"] == "A" :
                alt = l["AU_TUMOR"].split(",")[0]
            elif l["Alt"] == "C" :
                alt = l["CU_TUMOR"].split(",")[0]
            elif l["Alt"] == "G" :
                alt = l["GU_TUMOR"].split(",")[0]
            elif l["Alt"] == "T" :
                alt = l["TU_TUMOR"].split(",")[0]

            dp = float(alt) / (float(alt) + float(ref)) * 100
        elif somatic == "somatico" and "TAR" in l.keys() :
            dp = float(l["TIR_TUMOR"].split(",")[0]) / (float(l["TIR_TUMOR"].split(",")[0]) + float(l["TAR_TUMOR"].split(",")[0])) * 100

    except ZeroDivisionError :
        print "ERROR: La variant {} te cobertura 0".format([l["Chr"], l["Start"], l["End"], l["Ref"], l["Alt"]])

    return dp

def getMutationType(l, somatic) :
    ret = ""
    if somatic == "germinal" :
        ret = "GERMLINE"
    else :
        if "SOMATIC" in l.keys() :
            ret = l["SOMATIC"]
        else :
            ret = "SOMATIC"
    return ret

def maximMaf(dc) :
    maxim = -1
    for i in mafCols :
        if i in dc.keys() and dc[i] != 'NA' and dc[i] != '.' and float(dc[i]) > maxim:
            maxim = float(dc[i])
    if maxim == -1 :
        maxim = 'NA'

    return maxim

def resumPredictors(d) :
    deleterious = 0
    tolerated = 0
    unknown = 0
    if d["SIFT_pred"] == 'D' :
        deleterious += 1
    elif d["SIFT_pred"] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    if d["Polyphen2_HDIV_pred"] == 'D' or d["Polyphen2_HDIV_pred"] == 'P' :
        deleterious += 1
    elif d["Polyphen2_HDIV_pred"] == 'B' :
        tolerated += 1
    else :
        unknown += 1
    if d["Polyphen2_HVAR_pred"] == 'D' or d["Polyphen2_HVAR_pred"] == 'P' :
        deleterious += 1
    elif d["Polyphen2_HVAR_pred"] == 'B' :
        tolerated += 1
    else :
        unknown += 1
    if d["LRT_pred"] == 'D' :
        deleterious += 1
    elif d["LRT_pred"] == 'N' :
        tolerated += 1
    else :
        unknown += 1
    if d["MutationTaster_pred"] == 'A' or d["MutationTaster_pred"] == 'D' :
        deleterious += 1
    elif d["MutationTaster_pred"] == 'N' or d["MutationTaster_pred"] == 'P' :
        tolerated += 1
    else :
        unknown += 1
    if d["MutationAssessor_pred"] == 'H' or d["MutationAssessor_pred"] == 'M' :
        deleterious += 1
    elif d["MutationAssessor_pred"] == 'L' or d["MutationAssessor_pred"] == 'N' :
        tolerated += 1
    else :
        unknown += 1
    if d["FATHMM_pred"] == 'D' :
        deleterious += 1
    elif d["FATHMM_pred"] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    if d["PROVEAN_pred"] == 'D' :
        deleterious += 1
    elif d["PROVEAN_pred"] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    if d["MetaSVM_pred"] == 'D' :
        deleterious += 1
    elif d["MetaSVM_pred"] == 'T' :
        tolerated += 1
    else :
        unknown += 1
    if d["MetaLR_pred"] == 'D' :
        deleterious += 1
    elif d["MetaLR_pred"] == 'T' :
        tolerated += 1
    else :
        unknown += 1

    return "{}D, {}T, {}U".format(deleterious, tolerated, unknown)

def main(ruta) :
    # Leer el archivo multianno de ANNOVAR y guardar los datos en un diccionario
    convertirData(ruta)
    # Separar las variantes que han pasado todos los filtros de STrelka2
    # Guardar en un archivo de text todas las variantes (filtro0)
    # Separar las variantes exonicas (consecuencia) y las splicing en un diccionario aparte
    # Si el total de variantes es menor de 200 (numero arbitrario) re-anotar todas las variantes. En caso contrario, solo re-anotar las conseq
    # Agregar las columnas adicionales: population_max, predictor_summary, Strand_bias_score, GT, GQ, DP, RFD, ALD, VAF, MT, IGV, IGV_analisis, muestra
    # Guardar en un archivo de texto las variantes que pasan  (filtro0)
    # Guardar en un archivo de texto las variantes con consecuencia (filtro1)
    # Filtrar por MAF
    # Guardar en un archivo de texto las variantes con MAF>0.01 (filtro2)
    # Filtrar por VAF
    # Guardar en un archivo de texto las variantes con VAF <= 0.1 (filtro3)
    # Guardar en un archivo de texto las variantes con VAF > 0.1 (filtro4)

if __name__ == "__main__" :
    path = sys.argv[1]
    main(path)
"""
Columnas extra
l["population_max"] = maximMaf(l)
l["predictor_summary"] = resumPredictors(l)
l["Strand_bias_score"] = calcularStrandBias(l, somatic)
l["GT"] = getGenomeType(l, somatic)
l["GQ"] = getGenomeQuality(l, somatic)
l["DP"] = getRegionDepth(l, somatic)
l["RFD"] = getReferenceDepth(l, somatic)
l["ALD"] = getAlteratedDepth(l, somatic)
l["VAF"] = getVAF(l, somatic)
l["MT"] = getMutationType(l, somatic)
l["IGV"] = "http://localhost:60151/goto?locus={}:{}".format(l["Chr"], l["Start"])
l["IGV_analisis"] = ""
l["muestra"] = mostra
"""
