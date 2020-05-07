#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Recoge todos los datos del analisis y los guarda en un excel
"""

import os
import sys
import xlsxwriter

# Archivos que se usaran para recoger la informacion que se guardara en el excel final
qc = "../alnQC.txt"
cov = "../coverage.txt"
arxius = ["raw.reanno.tsv", "pass.reanno.tsv", "conseq.reanno.tsv", "highMAF.reanno.tsv". "lowVAF.reanno.tsv". "cand.reanno.tsv"]
stats = "variants.stats.txt"
# Orden de las columnas en que se colocaran en cada una de las pestanas del excel
orden = ["sample", "IGV_link", "Gene.refGene", "Chr", "Start", "End", "Ref", "Alt", "GT", "GQ", "GQX", "MQ", "Func.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene",
"Ref_depth", "Alt_depth", "DP", "DPF", "AD", "ADF", "ADR", "VAF", "population_max", "population_max_name", "predictor_summary", "Strand_bias_score", "SB",
"avsnp150", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "cosmic70",
"SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
"LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_pred", "FATHMM_score", "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred",
"VEST3_score", "MetaSVM_score", "MetaSVM_pred", "MetaLR_score", "MetaLR_pred", "M-CAP_score", "M-CAP_pred", "REVEL_score", "MutPred_score", "CADD_raw", "CADD_phred", "DANN_score",
"fathmm-MKL_coding_score", "fathmm-MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon_score", "integrated_fitCons_score", "integrated_confidence_value",
"GTEx_V6p_tissue", "GERP++_RS", "phyloP100way_vertebrate", "phyloP20way_mammalian", "phastCons100way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", "Interpro_domain", "GTEx_V6p_gene"]

# cabecera = ["Analysis", "Interpretation", "Sample", "IGV link", "IGV Analysis", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Genome type",
# "Mutation type", "Genome type", "Quality",
# "Depth (Normal sample, Tumor sample if applicable)", "Reference depth", "Alt depth", "VAF", "Function", "Transcript", "Exonic consequence", "AA change & transcript", "avsnp147", "1000g2015aug_all",
# "1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
# "esp6500siv2_all", "esp6500siv2_aa", "esp6500siv2_ea", "dbSNP_MAF", "CADD_1000g_all", "CADD_1000g_afr", "CADD_1000g_amr", "CADD_1000g_eas", "CADD_1000g_eur", "CADD_1000g_sas", "dbNSFP_1000g_all",
# "dbNSFP_1000g_afr", "dbNSFP_1000g_amr",	"dbNSFP_1000g_eas", "dbNSFP_1000g_eur", "dbNSFP_1000g_sas", "ExAC_ExAC_all", "ExAC_ExAC_afr", "ExAC_ExAC_amr", "ExAC_ExAC_eas", "ExAC_ExAC_fin", "ExAC_ExAC_nfe",
# "ExAC_ExAC_oth", "ExAC_ExAC_sas", "dbNSFP_ExAC_all", "dbNSFP_ExAC_afr", "dbNSFP_ExAC_amr", "dbNSFP_ExAC_eas", "dbNSFP_ExAC_fin", "dbNSFP_ExAC_nfe", "dbNSFP_ExAC_oth", "dbNSFP_ExAC_sas", "CADD_ESP6500_all",
# "CADD_ESP6500_aa", "CADD_ESP6500_ea", "dbNSFP_esp6500_all", "dbNSFP_esp6500_aa", "dbNSFP_esp6500_ea", "CLINSIG", "CLNACC", "CLNDBN", "CLNDSDB", "CLNDSDBID", "cosmic70", "SIFT_score", "SIFT_pred",
# "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score",
# "MutationAssessor_pred", "FATHMM_score", "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred", "VEST3_score", "CADD_raw", "CADD_phred", "DANN_score", "fathmm-MKL_coding_score", "fathmm-MKL_coding_pred",
# "MetaSVM_score", "MetaSVM_pred", "MetaLR_score", "MetaLR_pred", "integrated_fitCons_score", "integrated_confidence_value", "GERP++_RS", "phyloP7way_vertebrate", "phyloP20way_mammalian",
# "phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", "VCF FILTER", "INFO from VCF", "FORMAT from VCF", "FORMAT VALUES from VCF"]]

def convertirArchivo(path) :
    """
        Convierte un archivo .reanno.tsv en un diccionario
        Se asume que tiene la primera fila es la cabecera, la cual se usa para las claves del diccionario
    """
    cabecera = True
    datos = []
    with open(path, "r") as fi :
        for l in fi :
            if cabecera :
                claves = l.strip().split("\t")
                cabecera = False
            else :
                aux = l.strip().split("\t")
                it = 0
                temp = {}
                for c in claves :
                    temp[c] = aux[it]
                    it += 1
                datos.append(temp)
    return datos

def crearCabecera(hoja, libro) :
    """
    Escribe la cabecera de las pestanas de variantes
    """
    titulo = libro.add_format({'bold' : True,
        'align' : 'center',
        'border' : 1,
        'bg_color' :  '#B3E6FF',
        'font_size' : 13
    })
    # TODO: Escribir bien la cabecera
    cabecera = orden
    cols = 0
    nextLine = 0 # Filas donde se esta escribiendo la cabecera. Se devuelve a la funcion principal para no sobreescribir
    for n in cabecera :
        hoja.write(nextLine, cols, n, titulo)
        cols += 1
    return nextLine

def escribirVariantes(hoja, libro, cnt, empiezaEn) :
    """Escribir los datos del diccionario pasado por parametro en la pestana pasada por parametro
    hoja : hoja actual
    libro : libro excel donde se estan guardando los datos
    cnt : dict
        Diccionario con las variantes que se va a guardar en el excel
    empiezaEn : int
        Fila por la que se empezara a escribir el excel"""
    # Estilos para los predictores
    rojo = excel.add_format({"bg_color" : "#FF4D4D"})
    verde = excel.add_format({"bg_color" : "#43F906"})
    amarillo = excel.add_format({"bg_color" : "#FFFF00"})
    naranja = excel.add_format({"bg_color" : "#FF8000"})
    predictors = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred"]
    fila = empiezaEn
    columna = 0
    for dic in cnt :
        columna = 0
        for o in orden :
            if o in predictors :
                if dic[o] == 'D' : #Todos los predictores anotan una D como deleterea
                    hoja.write(fila, columna, dic[o], rojo)
                elif dic[o] == 'T' : #SIFT, Provean, MetaSVM, MetaLR y FATHMM anotan una T como tolerado
                    hoja.write(fila, columna, dic[o], verde)
                elif dic[o] == 'N' : #LRT, MutationTaster y MutationAssessor anotan una N como tolerado
                    hoja.write(fila, columna, dic[o], verde)
                elif (o == "Polyphen2_HDIV_pred" or o == "Polyphen2_HVAR_pred") and dic[o] == 'P' :
                    hoja.write(fila, columna, dic[o], naranja)
                elif (o == "Polyphen2_HDIV_pred" or o == "Polyphen2_HVAR_pred") and dic[o] == 'B' :
                    hoja.write(fila, columna, dic[o], verde)
                elif o == "MutationTaster_pred" and dic[o] == 'A' :
                    hoja.write(fila, columna, dic[o], rojo)
                elif o == "MutationTaster_pred" and dic[o] == 'P' :
                    hoja.write(fila, columna, dic[o], verde)
                elif o == "MutationAssessor_pred" and dic[o] == 'H' :
                    hoja.write(fila, columna, dic[o], rojo)
                elif o == "MutationAssessor_pred" and dic[o] == 'M' :
                    hoja.write(fila, columna, dic[o], naranja)
                elif o == "MutationAssessor_pred" and dic[o] == 'L' :
                    hoja.write(fila, columna, dic[o], verde)
                else :
                    hoja.write(fila, columna, dic[o])
            else :
                hoja.write(fila, columna, dic[o])
            columna += 1
        fila += 1
    return fila

def crearExcel(nom) :
    """Recoge la informacion de las variantes desde los archivos .reanno.tsv
    Convierte los datos de cada archivo en un diccionario
    Lo escribe en una pestana del archivo excel con el nombre pasado por parametro
    Recoge la informacion de calidad desde los archivos coverage.txt, alnQC.txt y variants.stats.txt
    Guarda toda esta informacion en una pestana del excel cuyo nombre se ha pasado por parametro
    """
    wb = xlsxwriter.Workbook("{}.xlsx".format(nom), {"strings_to_numbers" : True})
    full = wb.add_worksheet("Filtered_def") #Hoja vacio para colocar las variantes que pasan el filtro visual
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    #Comprobar si existen los archivos de variantes, leerlos y montarlos en el excel
    for a in arxius :
        if os.path.isfile(a) :
            aux = a.split(".")[0] # Recoger el nombre que tendra la pestana
            full = wb.add_worksheet(a)
            dc = convertirArchivo(a)
            filaActual = crearCabecera(full, wb)
            filaActual = escribirVariantes(full, wb, dc, filaActual)

if __name__ == "__main__" :
    # Comprobar que todos los archivos necesarios estan creados
    trobats = 0
    if os.path.isfile(qc) :
        trobats += 1
    else :
        print("WARNING: No encontrado el archivo con el control de calidad del alineamiento. Deberia estar en: {}".format(qc))
    if os.path.isfile(stats) :
        trobats += 1
    else :
        print("WARNING: No encontrado el archivo con las estadisticas de variantes. Buscado como: {}".format(stats))
    if os.path.isfile(cov) :
        trobats += 1
    else :
        print("WARNING: No encontrado el archivo con las estatidisticas de coverage. Buscado como: {}".format(cov))
    for a in arxius :
        if os.path.isfile(a) :
            trobats += 1
        else :
            print("WARNING: No encontrado el filtro de variantes {}".format(a))

    if trobats > 0 :
        if (len(sys.argv) > 1) :
            crearExcel(sys.argv[1])
        else :
            crearExcel("noName")
    else :
        print("ERROR: No se ha encontrado ninguno de los archivos necesarios. No se puede crear el excel")
