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
arxius = ["cand.reanno.tsv", "lowVAF.reanno.tsv", "highMAF.reanno.tsv", "conseq.reanno.tsv", "pass.reanno.tsv", "raw.reanno.tsv"]
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
    cabecera = ["Analysis", "Interpretation", "sample", "IGV_link", "IGV", "Gene.refGene", "Chr", "Start", "End", "Ref", "Alt", "GT", "GQ", "GQX", "MQ", "Func.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene",
    "Ref_depth", "Alt_depth", "DP", "DPF", "AD", "ADF", "ADR", "VAF", "population_max", "population_max_name", "predictor_summary", "Strand_bias_score", "SB",
    "avsnp150", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG", "cosmic70",
    "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
    "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score", "MutationAssessor_pred", "FATHMM_score", "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred",
    "VEST3_score", "MetaSVM_score", "MetaSVM_pred", "MetaLR_score", "MetaLR_pred", "M-CAP_score", "M-CAP_pred", "REVEL_score", "MutPred_score", "CADD_raw", "CADD_phred", "DANN_score",
    "fathmm-MKL_coding_score", "fathmm-MKL_coding_pred", "Eigen_coding_or_noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon_score", "integrated_fitCons_score", "integrated_confidence_value",
    "GTEx_V6p_tissue", "GERP++_RS", "phyloP100way_vertebrate", "phyloP20way_mammalian", "phastCons100way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", "Interpro_domain", "GTEx_V6p_gene"]
    cabecera = orden
    cols = 0
    nextLine = 0 # Filas donde se esta escribiendo la cabecera. Se devuelve a la funcion principal para no sobreescribir
    for n in cabecera :
        hoja.write(nextLine, cols, n, titulo)
        cols += 1
    nextLine += 1
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
    rojo = libro.add_format({"bg_color" : "#FF4D4D"})
    verde = libro.add_format({"bg_color" : "#43F906"})
    amarillo = libro.add_format({"bg_color" : "#FFFF00"})
    naranja = libro.add_format({"bg_color" : "#FF8000"})
    predictors = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred"]
    fila = empiezaEn
    columna = 2
    for dic in cnt :
        columna = 2
        for o in orden :
            if columna != 4 :
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

def ayudaPredictores(hoja, libro, fila) :
    """Escribe un cuadro con temas de ayuda interesantes para el usuario final a partir de la fila pasada como parametro. Los temas de ayuda son:
    Los predictores, los coeficientes de strand-bias, las calidades de la variante, dependiendo de si la muestra analizada es somatica o germinal, y la columna summary_predictors"""
    # Estilos
    titulo = libro.add_format({ 'bold' : True,
        'bg_color' : '#B3C6FF',
        'align' : 'center',
        'top' : 1,
        'left' : 1,
        'right' : 1,
        'font_size' : 13
    })
    medio = libro.add_format({
        'left' : 1,
        'right' : 1
    })
    bajo = libro.add_format({'bottom' : 1,
        'left' : 1,
        'right' : 1
    })

    #Ayuda para interpretar la prediccion de algunos de los predictores
    titol = "HELP for in silico predictors"
    sift = "SIFT\n\rD -> Deleterious\n\rT -> Benign"
    polyphen = "Polyphen\n\rD -> Probably damaging\n\rP -> Possibly damaging\n\rB -> Benign"
    lrt = "LRT\n\rD -> Deletereous\n\rN -> Neutral\n\rU -> Unknown"
    mtaster = "Mutation Taster\n\rA -> Disease causing automatic\n\rD -> Disease causing\n\rN -> Polymophism\n\rP -> Polymophism automatic"
    massessor = "Mutation Assessor\n\rH -> High (Deleterious)\n\rM -> Medium (Deleterious)\n\rL -> Low (Tolerated)\n\rN -> Neutral (Tolerated)"
    provean = "PROVEAN\n\rD -> Deleterious\n\rT -> Benign"
    fathmm = "FATHMM\n\rD -> Deleterious\n\rN -> Benign"
    msvm = "MetaSVM\n\rD -> Deleterious\n\rT -> Benign"
    mlr = "MetaLR\n\rD -> Deleterious\n\rT -> Benign"
    ot = "GERP++, PhyloP, SiPhy -> Higher scores are more deleterious"

    hoja.merge_range(fila, 2, fila, 8, titol, titulo)
    hoja.merge_range(fila+1, 2, fila+3, 8, sift, medio)
    hoja.merge_range(fila+4, 2, fila+7, 8, polyphen, medio)
    hoja.merge_range(fila+8, 2, fila+11, 8, lrt, medio)
    hoja.merge_range(fila+12, 2, fila+16, 8, mtaster, medio)
    hoja.merge_range(fila+17, 2, fila+21, 8, massessor, medio)
    hoja.merge_range(fila+22, 2, fila+24, 8, provean, medio)
    hoja.merge_range(fila+25, 2, fila+27, 8, fathmm, medio)
    hoja.merge_range(fila+28, 2, fila+30, 8, msvm, medio)
    hoja.merge_range(fila+31, 2, fila+33, 8, mlr, medio)
    hoja.merge_range(fila+34, 2, fila+34, 8, ot, bajo)

    #Ayuda para interpretar la columna summary_predictors
    titol = "HELP for 'predictor summary' column"
    body = "The column summarizes the prediction of SIFT,\n\rPolyphen2 HDIV, Polyphen2 HVAR, LRT, MutationTaster,\n\rMutationAssessor, FATHMM, PROVEAN, MetaSVM, and MetaLR.\n\n"
    body += "It enumerates the number of (D)eletereous, (T)olerated,\n\r and (U)nknown prediction\n\n"
    body += "So 2D, 7T, 1U\nmeans\n2 deleterious, 7 tolerated, and 1 unknown predictions."
    hoja.merge_range(fila, 10, fila, 16, titol, titulo)
    hoja.merge_range(fila+1, 10, fila+9, 16, body, bajo)

    #Ayuda para interpretar las calidades de la variante, dependiendo de si la muestra analizada es somatica o germinal
    titol = "HELP for quality column, depending on variant calling analysis"
    body = "Quality column represents the score for variant sites."

    hoja.merge_range(fila+12, 10, fila+12, 18, titol, titulo)
    hoja.merge_range(fila+13, 10, fila+14, 18, body, bajo)

    #Ayuda para los coeficientes de strand-bias
    titol = "HELP for Strand bias ratio"
    body = "Strand bias is calculated with the formula:\nreads_reference_forward + reads_alterated_forward - read_reference_reverse - reads_alterated_reverse"

    hoja.merge_range(fila+17, 10, fila+17, 19, titol, titulo)
    hoja.merge_range(fila+18, 10, fila+19, 19, body, bajo)

def escribirEstadisticas(hoja, libro) :
    if os.path.isfile(qc) :
        hoja.write(0, 0, "Calidad del alineamiento")
    if os.path.isfile(cov) :
        hoja.write(0, 4, "Datos de coverage")
    if os.path.isfile(stats) :
        hoja.write(0, 8, "Estadisticas de variantes")

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
            full = wb.add_worksheet(aux.capitalize())
            dc = convertirArchivo(a)
            filaActual = crearCabecera(full, wb)
            filaActual = escribirVariantes(full, wb, dc, filaActual)
            ayudaPredictores(full, wb, filaActual+2)

    full = wb.add_worksheet("QC_stats")
    escribirEstadisticas(full, wb)

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
