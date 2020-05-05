# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
MAIN: Re-anotar y filtrar el archivo de salida de table_annovar
"""
import sys
import getWebInfo as gw

#Todas las claves con la informacion que se guardara en los archivos de texto
allkeys = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "1000g2015aug_all",
"1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH",
"ExAC_SAS", "gnomad_exome_AF", "gnomad_exome_AF_popmax", "gnomad_exome_AF_male", "gnomad_exome_AF_female", "gnomad_exome_AF_raw", "gnomad_exome_AF_afr", "gnomad_exome_AF_sas",
"gnomad_exome_AF_amr", "gnomad_exome_AF_eas", "gnomad_exome_AF_nfe", "gnomad_exome_AF_fin", "gnomad_exome_AF_asj", "gnomad_exome_AF_oth", "gnomad_exome_non_topmed_AF_popmax",
"gnomad_exome_non_neuro_AF_popmax", "gnomad_exome_non_cancer_AF_popmax", "gnomad_exome_controls_AF_popmax", "gnomad_genome_AF", "gnomad_genome_AF_popmax", "gnomad_genome_AF_male",
"gnomad_genome_AF_female", "gnomad_genome_AF_raw", "gnomad_genome_AF_afr", "gnomad_genome_AF_sas", "gnomad_genome_AF_amr", "gnomad_genome_AF_eas", "gnomad_genome_AF_nfe",
"gnomad_genome_AF_fin", "gnomad_genome_AF_asj", "gnomad_genome_AF_oth", "gnomad_genome_non_topmed_AF_popmax", "gnomad_genome_non_neuro_AF_popmax", "gnomad_genome_non_cancer_AF_popmax",
"gnomad_genome_controls_AF_popmax", "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "cosmic70", "SIFT_score",
"SIFT_converted_rankscore", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore", "Polyphen2_HVAR_pred",
"LRT_score", "LRT_converted_rankscore", "LRT_pred", "MutationTaster_score", "MutationTaster_converted_rankscore", "MutationTaster_pred", "MutationAssessor_score",
"MutationAssessor_score_rankscore", "MutationAssessor_pred", "FATHMM_score", "FATHMM_converted_rankscore", "FATHMM_pred", "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred",
"VEST3_score", "VEST3_rankscore", "MetaSVM_score", "MetaSVM_rankscore", "MetaSVM_pred", "MetaLR_score", "MetaLR_rankscore", "MetaLR_pred", "M-CAP_score", "M-CAP_rankscore", "M-CAP_pred",
"REVEL_score", "REVEL_rankscore", "MutPred_score", "MutPred_rankscore", "CADD_raw", "CADD_raw_rankscore", "CADD_phred", "DANN_score", "DANN_rankscore", "fathmm-MKL_coding_score",
"fathmm-MKL_coding_rankscore", "fathmm-MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon_score", "GenoCanyon_score_rankscore", "integrated_fitCons_score",
"integrated_fitCons_score_rankscore", "integrated_confidence_value", "GERP++_RS", "GERP++_RS_rankscore", "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore",
"phyloP20way_mammalian", "phyloP20way_mammalian_rankscore", "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore",  "phastCons20way_mammalian",
"phastCons20way_mammalian_rankscore", "SiPhy_29way_logOdds", "SiPhy_29way_logOdds_rankscore", "Interpro_domain", "GTEx_V6p_gene", "GTEx_V6p_tissue", "CHROM", "POS", "ID", "REFERENCE",
"ALTERATED", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", "END", "BLOCKAVG", "SNVHPOL", "CIGAR", "RU", "REFPREP", "IDREP", "MQ", "GT", "GQ", "GQX", "DP", "DPF", "MIN_DP", "AD",
"ADF", "ADR", "FT", "DPI", "PL", "PS", "SB", "population_max", "predictor_summary", "Strand_bias_score", "Ref_depth", "Alt_depth", "VAF", "IGV_link", "sample", "CADD_1000g_all",
"CADD_1000g_afr", "CADD_1000g_amr", "CADD_1000g_eur", "CADD_1000g_eas", "CADD_1000g_sas", "dbNSFP_1000g_all", "dbNSFP_1000g_afr", "dbNSFP_1000g_amr", "dbNSFP_1000g_eur", "dbNSFP_1000g_eas",
"dbNSFP_1000g_sas", "CADD_ESP6500_all", "CADD_ESP6500_ea", "CADD_ESP6500_aa","dbNSFP_esp6500_all", "dbNSFP_esp6500_ea", "dbNSFP_esp6500_aa", "ExAC_ExAC_all", "ExAC_ExAC_afr","ExAC_ExAC_amr",
"ExAC_ExAC_eas", "ExAC_ExAC_fin", "ExAC_ExAC_nfe", "ExAC_ExAC_oth", "ExAC_ExAC_sas", "dbNSFP_ExAC_all", "dbNSFP_ExAC_afr", "dbNSFP_ExAC_amr", "dbNSFP_ExAC_eas", "dbNSFP_ExAC_fin",
"dbNSFP_ExAC_nfe", "dbNSFP_ExAC_oth", "dbNSFP_ExAC_sas", "gNOMAD_Exome_all", "gNOMAD_Exome_afr", "gNOMAD_Exome_amr", "gNOMAD_Exome_asj", "gNOMAD_Exome_eas", "gNOMAD_Exome_fin",
"gNOMAD_Exome_nfe", "gNOMAD_Exome_oth", "gNOMAD_Exome_popmax", "gNOMAD_Exome_raw", "gNOMAD_Exome_sas", "gNOMAD_Genome_all", "gNOMAD_Genome_afr", "gNOMAD_Genome_amr", "gNOMAD_Genome_asj",
"gNOMAD_Genome_eas", "gNOMAD_Genome_fin", "gNOMAD_Genome_nfe", "gNOMAD_Genome_oth", "gNOMAD_Genome_popmax", "gNOMAD_Genome_raw", "dbSNP_MAF"]

# Cuando un dato no este disponible se guarda el siguinte caracter
vacio = "--"

#Columnas que contienen informacion de MAF
mafCols = ["1000g2015aug_all", "1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN",
"ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "gnomad_exome_AF", "gnomad_exome_AF_popmax", "gnomad_exome_AF_male", "gnomad_exome_AF_female", "gnomad_exome_AF_raw", "gnomad_exome_AF_afr",
"gnomad_exome_AF_sas", "gnomad_exome_AF_amr", "gnomad_exome_AF_eas", "gnomad_exome_AF_nfe", "gnomad_exome_AF_fin", "gnomad_exome_AF_asj", "gnomad_exome_AF_oth", "gnomad_exome_non_topmed_AF_popmax",
"gnomad_exome_non_neuro_AF_popmax", "gnomad_exome_non_cancer_AF_popmax", "gnomad_exome_controls_AF_popmax", "gnomad_genome_AF", "gnomad_genome_AF_popmax", "gnomad_genome_AF_male",
"gnomad_genome_AF_female", "gnomad_genome_AF_raw", "gnomad_genome_AF_afr", "gnomad_genome_AF_sas", "gnomad_genome_AF_amr", "gnomad_genome_AF_eas", "gnomad_genome_AF_nfe",
"gnomad_genome_AF_fin", "gnomad_genome_AF_asj", "gnomad_genome_AF_oth", "gnomad_genome_non_topmed_AF_popmax", "gnomad_genome_non_neuro_AF_popmax", "gnomad_genome_non_cancer_AF_popmax",
"gnomad_genome_controls_AF_popmax", "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "CADD_1000g_all", "CADD_1000g_afr", "CADD_1000g_amr", "CADD_1000g_eur", "CADD_1000g_eas",
"CADD_1000g_sas", "dbNSFP_1000g_all", "dbNSFP_1000g_afr", "dbNSFP_1000g_amr", "dbNSFP_1000g_eur", "dbNSFP_1000g_eas", "dbNSFP_1000g_sas", "CADD_ESP6500_all", "CADD_ESP6500_ea", "CADD_ESP6500_aa",
"dbNSFP_esp6500_all", "dbNSFP_esp6500_ea", "dbNSFP_esp6500_aa", "ExAC_ExAC_all", "ExAC_ExAC_afr","ExAC_ExAC_amr", "ExAC_ExAC_eas", "ExAC_ExAC_fin", "ExAC_ExAC_nfe", "ExAC_ExAC_oth",
"ExAC_ExAC_sas", "dbNSFP_ExAC_all", "dbNSFP_ExAC_afr", "dbNSFP_ExAC_amr", "dbNSFP_ExAC_eas", "dbNSFP_ExAC_fin", "dbNSFP_ExAC_nfe", "dbNSFP_ExAC_oth", "dbNSFP_ExAC_sas", "gNOMAD_Exome_all",
"gNOMAD_Exome_afr", "gNOMAD_Exome_amr", "gNOMAD_Exome_asj", "gNOMAD_Exome_eas", "gNOMAD_Exome_fin", "gNOMAD_Exome_nfe", "gNOMAD_Exome_oth", "gNOMAD_Exome_popmax", "gNOMAD_Exome_raw",
"gNOMAD_Exome_sas", "gNOMAD_Genome_all", "gNOMAD_Genome_afr", "gNOMAD_Genome_amr", "gNOMAD_Genome_asj", "gNOMAD_Genome_eas", "gNOMAD_Genome_fin", "gNOMAD_Genome_nfe", "gNOMAD_Genome_oth",
"gNOMAD_Genome_popmax", "gNOMAD_Genome_raw", "dbSNP_MAF"]

def addFORMAT(dic) :
    """
    Separar los datos de las columnas INFO y FORMAT del vcf en columnas. Se devuelven los datos separados en un diccionario. Se comprueba si alguna clave esta ya presente en el diccionario
    """
    claves = dic["FORMAT"].split(":")
    valores = dic["SAMPLE"].split(":")
    minidic = {}
    i = 0
    #Separar los datos de la columna FORMAT en distintas columnas
    for c in claves :
        if c in dic.keys() :
            print("WARNING: Clau {} repetida en {}".format(c, dic.keys()))
            sys.exit()
        else :
            minidic[c] = valores[i]
        i += 1

    #Separar los datos de la columna INFO en distintas columnas
    datos = dic["INFO"].split(";")
    for d in datos :
        aux = d.split("=")
        aux[0] = "{}_INFO".format(aux[0])
        if aux[0] in dic.keys() or aux[0] in minidic.keys():
            print("WARNING: Clau {} repetida al muntar INFO".format(aux[0]))
            sys.exit()
        if len(aux) == 1 :
            minidic[aux[0]] = aux[0]
        else :
            minidic[aux[0]] = aux[1]
    return minidic

def convertirData(path) :
    """Abrir un archivo con formato table_annovar y guardar todo el contenido en un diccionario"""
    cabecera = True
    temp = {}
    # Columnas dentro del table_annovar
    claves = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "1000g2015aug_all",
    "1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH",
    "ExAC_SAS", "gnomad_exome_AF", "gnomad_exome_AF_popmax", "gnomad_exome_AF_male", "gnomad_exome_AF_female", "gnomad_exome_AF_raw", "gnomad_exome_AF_afr", "gnomad_exome_AF_sas",
    "gnomad_exome_AF_amr", "gnomad_exome_AF_eas", "gnomad_exome_AF_nfe", "gnomad_exome_AF_fin", "gnomad_exome_AF_asj", "gnomad_exome_AF_oth", "gnomad_exome_non_topmed_AF_popmax",
    "gnomad_exome_non_neuro_AF_popmax", "gnomad_exome_non_cancer_AF_popmax", "gnomad_exome_controls_AF_popmax", "gnomad_genome_AF", "gnomad_genome_AF_popmax", "gnomad_genome_AF_male",
    "gnomad_genome_AF_female", "gnomad_genome_AF_raw", "gnomad_genome_AF_afr", "gnomad_genome_AF_sas", "gnomad_genome_AF_amr", "gnomad_genome_AF_eas", "gnomad_genome_AF_nfe",
    "gnomad_genome_AF_fin", "gnomad_genome_AF_asj", "gnomad_genome_AF_oth", "gnomad_genome_non_topmed_AF_popmax", "gnomad_genome_non_neuro_AF_popmax", "gnomad_genome_non_cancer_AF_popmax",
    "gnomad_genome_controls_AF_popmax", "esp6500siv2_all", "esp6500siv2_ea", "esp6500siv2_aa", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "cosmic70", "SIFT_score",
    "SIFT_converted_rankscore", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore", "Polyphen2_HVAR_pred",
    "LRT_score", "LRT_converted_rankscore", "LRT_pred", "MutationTaster_score", "MutationTaster_converted_rankscore", "MutationTaster_pred", "MutationAssessor_score",
    "MutationAssessor_score_rankscore", "MutationAssessor_pred", "FATHMM_score", "FATHMM_converted_rankscore", "FATHMM_pred", "PROVEAN_score", "PROVEAN_converted_rankscore", "PROVEAN_pred",
    "VEST3_score", "VEST3_rankscore", "MetaSVM_score", "MetaSVM_rankscore", "MetaSVM_pred", "MetaLR_score", "MetaLR_rankscore", "MetaLR_pred", "M-CAP_score", "M-CAP_rankscore", "M-CAP_pred",
    "REVEL_score", "REVEL_rankscore", "MutPred_score", "MutPred_rankscore", "CADD_raw", "CADD_raw_rankscore", "CADD_phred", "DANN_score", "DANN_rankscore", "fathmm-MKL_coding_score",
    "fathmm-MKL_coding_rankscore", "fathmm-MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon_score", "GenoCanyon_score_rankscore", "integrated_fitCons_score",
    "integrated_fitCons_score_rankscore", "integrated_confidence_value", "GERP++_RS", "GERP++_RS_rankscore", "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore",
    "phyloP20way_mammalian", "phyloP20way_mammalian_rankscore", "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore",  "phastCons20way_mammalian",
    "phastCons20way_mammalian_rankscore", "SiPhy_29way_logOdds", "SiPhy_29way_logOdds_rankscore", "Interpro_domain", "GTEx_V6p_gene", "GTEx_V6p_tissue", "CHROM", "POS", "ID", "REFERENCE",
    "ALTERATED", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    with open(path, "r") as fi :
        dc = []
        for l in fi :
            if cabecera :
                cabecera = False
            else :
                aux = l.split("\t")
                i = 0
                for c in claves :
                    temp[c] = aux[i]
                    i += 1
                aux2 = addFORMAT(temp)
                temp.update(aux2)
                dc.append(temp)
                temp = {}
    return dc

def filtroPAS(dc) :
    """Seleccionar solo las variantes que tengan como filtro PASS segun Strelka2"""
    nuevo = []
    for d in dc :
        if d["FILTER"] == "PASS" :
            nuevo.append(d)
    return nuevo

def filtroConseq(dc) :
    """Selecciona aquellas variantes que son exonicas (excluye sinonimas) y splicing"""
    nuevo = []
    for d in dc :
        if d["Func.refGene"] == "splicing" :
            nuevo.append(d)
        elif d["Func.refGene"] == "exonic" and d["ExonicFunc.refGene"] != "synonymous SNV" :
            nuevo.append(d)
    return nuevo

def filtrarMAF(dc) :
    """Clasificar las variantes por su Minor Allele Frequency"""
    baja = []
    alta = []
    for l in dc :
        if l["population_max"] != "NA" and l["population_max"] >= 0.01 :
            alta.append(l)
        else :
            baja.append(l)
    return alta, baja

def filtrarVAF(dc) :
    """Clasificar las variantes por su Variant Allele Frequency"""
    alta = []
    baja = []
    for l in dc :
        if l["VAF"] > 10 :
            alta.append(l)
        else :
            baja.append(l)
    return alta, baja

def addWebInfo(dc) :
    """Agregar la informacion de myvariant.info a aquellas variantes que sean exonicas"""
    for d in dc :
        aux = gw.getMyVariant(d["Chr"], d["Start"], d["Ref"], d["Alt"])
        d.update(aux)
    return dc

def guardarTabla(dc, prefijo) :
    """Guardar el diccionario pasado por parametro en un archivo de texto con formato de tabla. El prefijo es el nombre que tendra el archivo"""
    sufijo = "reanno.tsv"
    filename = "{}.{}".format(prefijo, sufijo)
    with open(filename, "w") as fi :
        fi.write("\t".join(allkeys))
        fi.write("\n")
        for d in dc :
            for o in allkeys :
                if o in d.keys() :
                    fi.write("{}\t".format(d[o]))
                else :
                    fi.write("{}\t".format(vacio))
            fi.write("\n")
    print("INFO: Guardado archivo {}".format(filename))


def calcularStrandBias(l) :
    """Calcular el numero total de reads en forward - numero total de reads en reverse"""
    coef = "NP"
    auxF = l["ADF"].split(",")
    auxR = l["ADR"].split(",")
    forward = 0
    reverse = 0
    for a in auxF :
        forward += int(a)
    for a in auxR :
        reverse += int(a)
    coef = forward - reverse

    return str(coef)

def maximMaf(dc) :
    maxim = -1
    poblacio = ""
    for i in mafCols :
        if i in dc.keys() and dc[i] != 'NA' and dc[i] != '.' and float(dc[i]) > maxim:
            maxim = float(dc[i])
            maxim = i
    if maxim == -1 :
        maxim = 'NA'
    else :
        print("La poblacio escollida es: {}".format(poblacio))
    return maxim, poblacio

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

def main(ruta, samplename = "noName") :
    webInfo = False
    # Leer el archivo multianno de ANNOVAR y guardar los datos en un diccionario
    todas = convertirData(ruta)
    # Separar las variantes que han pasado todos los filtros de STrelka2
    pas = filtroPAS(todas)
    if len(pas) <= 200 :
        pas = addWebInfo(pas)
        webInfo = True
    # Agregar las columnas adicionales: "population_max", "predictor_summary". "Strand_bias_score", "Ref_depth", "Alt_depth", "VAF", "IGV_link", "sample"
    for p in pas :
        p["population_max"], p["population_max_name"] = maximMaf(p)
        p["predictor_summary"] = resumPredictors(p)
        p["Strand_bias_score"] = calcularStrandBias(p)
        p["Ref_depth"] = p["AD"].split(",")[0]
        p["Alt_depth"] = p["AD"].split(",")[1]
        if "DP" in p.keys() :
            p["VAF"] = float(p["Alt_depth"])/float((p["DP"])) * 100
        else :
            tmp = p["AD"].split(",")
            dp = 0
            for t in tmp :
                dp += int(t)
            p["VAF"] = float(p["Alt_depth"])/dp * 100
        p["IGV_link"] = "http://localhost:60151/goto?locus={}:{}".format(p["Chr"], p["Start"])
        p["sample"] = samplename

    # Separar las variantes exonicas (consecuencia) y las splicing en un diccionario aparte
    conseq = filtroConseq(pas)
    # Si el total de variantes es menor de 200 (numero arbitrario) re-anotar todas las variantes. En caso contrario, solo re-anotar las conseq
    if not webInfo :
        conseq = addWebInfo(conseq)
    # Filtrar por MAF
    mafAlta, mafBaja = filtrarMAF(conseq)
    # Filtrar por VAF
    vafAlta, vafBaja = filtrarVAF(mafBaja)

    # Guardar los filtros en archivos de texto
    guardarTabla(todas, "raw") # Todas las variantes recogidas en el vcf
    guardarTabla(pas, "pass") # Variantes que Strelka2 considera que han pasado todos sus filtros
    guardarTabla(conseq, "conseq") # Variantes en regiones de splicing o exonicas (exluyendo sinonimas)
    guardarTabla(mafAlta, "highMAF") # Variantes con una MAF>=0.01 en alguna columna de base de datos poblacional
    guardarTabla(vafBaja, "lowVAF") # Variantes con una frecuencia alelica menor de 0.1
    guardarTabla(vafAlta, "cand") # Variantes que han pasado todos los filtros anteriores

    # Mostrar por pantalla estadisticas basicas
    print("Totes les variants {}".format(len(todas))) # Total de variantes reportadas por Strelka2
    print("Variants amb filtres=PASS {}".format(len(pas)))
    print("Variants conseq: {}".format(len(conseq)))
    print("MAF >= 0.1: {}".format(len(mafAlta)))
    print("MAF < 0.1: {}".format(len(mafBaja)))
    print("Baixa VAF: {}".format(len(vafBaja)))
    print("Variants candidates: {}".format(len(vafAlta)))

    # Agregar estadisticas de los tipos de variantes recogidos en el panel a la pestaña QC
    with open("variants.stats.tsv", "w") as fi :
        fi.write("{")
        fi.write("Totales : {},".format(len(todas)))
        fi.write("Filtro=PASS : {},".format(len(pas)))
        fi.write("Conseq : {},".format(len(conseq)))
        fi.write("MAF_alta : {},".format(len(mafAlta)))
        fi.write("VAF_baja : {},".format(len(vafBaja)))
        fi.write("Candidatas : {}".format(len(vafAlta)))
        fi.write("}")

    with open("variants.stats.tsv", "r") as fi :
        aux = fi.read()

    d = eval(aux)
    print(d)
    print(d["Totales"])

if __name__ == "__main__" :
    path = sys.argv[1]
    sample = sys.argv[1]
    main(path, sample)
