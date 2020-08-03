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
covGens = "../coverage/coverageGeneStats.tsv"
arxius = ["cand.reanno.tsv", "lowVAF.reanno.tsv", "highMAF.reanno.tsv", "conseq.reanno.tsv", "raw.reanno.tsv"]
stats = "variants.stats.txt"
# Orden de las columnas en que se colocaran en cada una de las pestanas del excel
orden = ["sample", "IGV_link", "Gene.refGene", "Chr", "Start", "End", "Ref", "Alt", "GT", "GQ", "MQ", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene", "GeneDetail.refGene",
"Ref_depth", "Alt_depth", "DP", "AD", "ADF", "ADR", "VAF", "FILTER", "population_max", "population_max_name", "gnomad_exome_AF_popmax", "gnomad_exome_non_topmed_AF_popmax",
"gnomad_genome_AF_popmax", "gnomad_genome_non_topmed_AF_popmax",
"predictor_summary", "Strand_bias_score", "SB",
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

    cabecera = ["Analysis", "Interpretation", "Sample", "IGV_link", "IGV", "Gene", "Chromosome", "Start", "End", "Ref", "Alt", "Genotype", "Genome Quality", "Mapping Quality", "Mutation type",
    "Exonic mutation type", "Amino acid change", "Transcript",
    "Reference depth", "Alterated depth", "Position depth", "Allelic depths", "Forward allelic depths", "Reverse allelic depths", "VAF", "Variant caller filter", "Population max MAF",
    "Population reported max MAF", "gNOMAD Exome popmax", "gNOMAD Exome nonTOPMed popmax", "gNOMAD Genome popmax", "gNOMAD Genome nonTOPMed popmax",
    "Predictor summary", "SMD strand bias score", "Variant caller strand bias score",
    "DBSNP", "ClinVar CLNALLELEID", "ClinVar CLNDN", "ClinVar CLNDISDB", "ClinVar CLNREVSTAT", "Clinvar CLNSIG", "COSMIC",
    "SIFT score", "SIFT pred", "Polyphen2 HDIV score", "Polyphen2 HDIV pred", "Polyphen2 HVAR score", "Polyphen2 HVAR pred",
    "LRT score", "LRT pred", "MutationTaster score", "MutationTaster pred", "MutationAssessor score", "MutationAssessor pred", "FATHMM score", "FATHMM pred", "PROVEAN score", "PROVEAN pred",
    "VEST3 score", "MetaSVM score", "MetaSVM pred", "MetaLR score", "MetaLR pred", "M-CAP score", "M-CAP pred", "REVEL score", "MutPred score", "CADD raw", "CADD phred", "DANN score",
    "fathmm-MKL coding score", "fathmm-MKL coding pred", "Eigen coding or noncoding", "Eigen-raw", "Eigen-PC-raw", "GenoCanyon score", "integrated fitCons score", "integrated confidence value",
    "GTEx V6p tissue", "GERP++ RS", "phyloP100way vertebrate", "phyloP20way mammalian", "phastCons100way vertebrate", "phastCons20way mammalian", "SiPhy 29way logOdds", "Interpro domain", "GTEx V6p gene"]

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
            if columna == 4 :
                columna += 1
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
            elif o == "Chr" and not dic[o].startswith("chr") :
                hoja.write(fila, columna, "chr{}".format(dic[o]))
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

    hoja.merge_range(fila, 1, fila, 8, titol, titulo)
    hoja.merge_range(fila+1, 1, fila+3, 8, sift, medio)
    hoja.merge_range(fila+4, 1, fila+7, 8, polyphen, medio)
    hoja.merge_range(fila+8, 1, fila+11, 8, lrt, medio)
    hoja.merge_range(fila+12, 1, fila+16, 8, mtaster, medio)
    hoja.merge_range(fila+17, 1, fila+21, 8, massessor, medio)
    hoja.merge_range(fila+22, 1, fila+24, 8, provean, medio)
    hoja.merge_range(fila+25, 1, fila+27, 8, fathmm, medio)
    hoja.merge_range(fila+28, 1, fila+30, 8, msvm, medio)
    hoja.merge_range(fila+31, 1, fila+33, 8, mlr, medio)
    hoja.merge_range(fila+34, 1, fila+34, 8, ot, bajo)

    #Ayuda para interpretar la columna summary_predictors
    titol = "HELP for 'predictor summary' column"
    body = "The column summarizes the prediction of SIFT,\n\rPolyphen2 HDIV, Polyphen2 HVAR, LRT, MutationTaster,\n\rMutationAssessor, FATHMM, PROVEAN, MetaSVM, and MetaLR.\n\n"
    body += "It enumerates the number of (D)eletereous, (T)olerated,\n\r and (U)nknown prediction\n\n"
    body += "So 2D, 7T, 1U\nmeans\n2 deleterious, 7 tolerated, and 1 unknown predictions."
    hoja.merge_range(fila, 10, fila, 16, titol, titulo)
    hoja.merge_range(fila+1, 10, fila+9, 16, body, bajo)

    #Ayuda para los coeficientes de strand-bias
    titol = "HELP for SMD Strand bias score"
    body = "Strand bias is calculated with the formula:\nreads_reference_forward + reads_alterated_forward - read_reference_reverse - reads_alterated_reverse"

    hoja.merge_range(fila+12, 10, fila+12, 19, titol, titulo)
    hoja.merge_range(fila+13, 10, fila+14, 19, body, bajo)

def escribirEstadisticas(hoja, libro) :
    # Estilos
    titulo = libro.add_format({ 'bold' : True,
        'bg_color' : '#B3FFC6',
        'align' : 'center',
        'top' : 1,
        'left' : 1,
        'right' : 1,
        'font_size' : 13
    })
    izquierda = libro.add_format({
        'left' : 1,
        'bold' : True
    })
    centro = libro.add_format({
        'bold' : True
    })
    derecha = libro.add_format({
        'right' : 1
    })
    derechaTit = libro.add_format({
        'right' : 1,
        'bold' : True
    })
    bajoIzq = libro.add_format({'bottom' : 1,
        'left' : 1,
        'bold' : True
    })
    bajoCtr = libro.add_format({'bottom' : 1
    })
    bajoDrch = libro.add_format({'bottom' : 1,
        'right' : 1
    })
    if os.path.isfile(qc) :
        with open(qc, "r") as fi :
            aux = fi.read()
        qua = eval(aux)
        hoja.merge_range(0, 0, 0, 1, "Alignment quality", titulo)
        hoja.write(1, 0, "FASTQ reads", izquierda)
        hoja.write(1, 1, "{}".format(qua["FASTQ"]), derecha)
        hoja.write(2, 0, "Aligned reads", izquierda)
        hoja.write(2, 1, "{} ({:.2f} %)".format(qua["BAM"], 100*float(qua["BAM"])/float(qua["FASTQ"])), derecha)
        hoja.write(3, 0, "ON target", izquierda)
        hoja.write(3, 1, "{} ({:.2f} %)".format(qua["ON"], 100*float(qua["ON"])/float(qua["BAM"])), derecha)
        if "DUPS" in qua.keys() :
            hoja.write(4, 0, "OFF target", izquierda)
            hoja.write(4, 1, "{} ({:.2f} %)".format(qua["OFF"], 100*float(qua["OFF"])/float(qua["BAM"])), derecha)
            hoja.write(5, 0, "Duplicates", bajoIzq)
            hoja.write(5, 1, "{}".format(qua["DUPS"]), bajoDrch)
        else :
            hoja.write(4, 0, "OFF target", bajoIzq)
            hoja.write(4, 1, "{} ({:.2f} %)".format(qua["OFF"], 100*float(qua["OFF"])/float(qua["BAM"])), bajoDrch)
        # Incrustar las imagenes de calidad reportadas por FastQC
        if os.path.isdir("../fastqc") :
            basequal = []
            seqcontent = []
            for root, dirs, files in os.walk("../fastqc", topdown = False) :
                for n in files :
                    if n == "per_base_sequence_content.png" :
                        seqcontent.append(os.path.join(root, n))
                    if n == "per_base_quality.png" :
                        basequal.append(os.path.join(root, n))
            fila = 6
            for i in basequal :
                hoja.insert_image(fila, 0, i, {"x_scale" : 0.2, "y_scale" : 0.2})
                fila += 6

            for i in seqcontent :
                hoja.insert_image(fila, 0, i, {"x_scale" : 0.2, "y_scale" : 0.2})
                fila += 6


    if os.path.isfile(cov) :
        with open(cov, "r") as fi:
            aux = fi.read()
        cv = eval(aux)
        hoja.merge_range(0, 4, 0, 5, "Coverage statistics", titulo)
        hoja.write(1, 4, "Min", izquierda)
        hoja.write(1, 5, cv["minimo"], derecha)
        hoja.write(2, 4, "Max", izquierda)
        hoja.write(2, 5, cv["maximo"], derecha)
        hoja.write(3, 4, "Average", izquierda)
        hoja.write(3, 5, cv["media"], derecha)
        hoja.write(4, 4, "Median", izquierda)
        hoja.write(4, 5, cv["mediana"], derecha)
        hoja.write(5, 4, "% bases with 0 cov", izquierda)
        hoja.write(5, 5, cv["bases0"], derecha)
        hoja.write(6, 4, "% bases with 30 cov", izquierda)
        hoja.write(6, 5, cv["bases30"], derecha)
        hoja.write(7, 4, "% bases with 100 cov", izquierda)
        hoja.write(7, 5, cv["bases100"], derecha)
        hoja.write(8, 4, "% bases with 500 cov", izquierda)
        hoja.write(8, 5, cv["bases500"], derecha)
        hoja.write(9, 4, "% bases with 1000 cov", bajoIzq)
        hoja.write(9, 5, cv["bases1000"], bajoDrch)
        if os.path.isfile("../coverage/coverage.png") : #Adjuntar el grafico de coverage general. Se reduce un 75%
            hoja.insert_image(12, 4, "../coverage/coverage.png", {"x_scale" : 0.25, "y_scale" : 0.25})

    if os.path.isfile(covGens) :
        hoja.merge_range(0, 8, 0, 12, "Coverage per gene", titulo)
        header = True
        fila = 1
        with open(covGens, "r") as fi :
            for l in fi :
                aux = l.strip("\n").split("\t")
                columna = 8
                if header :
                    hoja.write(fila, 8, aux[0].strip("\"").capitalize(), izquierda) # Columna nombre del gen
                    hoja.write(fila, 9, aux[1].strip("\"").capitalize(), centro) # Columna minim
                    hoja.write(fila, 10, aux[2].strip("\"").capitalize(), centro) # Columna maxim
                    hoja.write(fila, 11, aux[3].strip("\"").capitalize(), centro) # Columna mitjana
                    hoja.write(fila, 12, aux[4].strip("\"").capitalize(), derechaTit) # Columna mediana
                    header = False
                else :
                    hoja.write(fila, 8, aux[0].strip("\""), izquierda)
                    hoja.write(fila, 9, aux[1])
                    hoja.write(fila, 10, aux[2])
                    hoja.write(fila, 11, "{:.2f}".format(float(aux[3])))
                    hoja.write(fila, 12, aux[4], derecha)
                fila += 1
            hoja.write(fila-1, 8, aux[0].strip("\""), bajoIzq)
            hoja.write(fila-1, 9, aux[1], bajoCtr)
            hoja.write(fila-1, 10, aux[2], bajoCtr)
            hoja.write(fila-1, 11, "{:.2f}".format(float(aux[3])), bajoCtr)
            hoja.write(fila-1, 12, aux[4], bajoDrch)


    if os.path.isfile(stats) :
        with open(stats, "r") as fi :
            aux = fi.read()
        var = eval(aux)
        hoja.merge_range(0, 14, 0, 15, "Variant number per filter", titulo)
        hoja.write(1, 14, "RAW", izquierda)
        hoja.write(1, 15, var["Totales"], derecha)
        del var["Totales"]
        hoja.write(2, 14, "Consequence", izquierda)
        hoja.write(2, 15, var["Conseq"], derecha)
        del var["Conseq"]
        hoja.write(3, 14, "High MAF", izquierda)
        hoja.write(3, 15, var["MAF_alta"], derecha)
        del var["MAF_alta"]
        hoja.write(4, 14, "Low VAF", izquierda)
        hoja.write(4, 15, var["VAF_baja"], derecha)
        del var["VAF_baja"]
        hoja.write(5, 14, "Candidates", bajoIzq)
        hoja.write(5, 15, var["Candidatas"], bajoDrch)
        del var["Candidatas"]
        hoja.merge_range(0, 17, 0, 18, "Raw variant distribution", titulo)
        fila = 1
        for k, v in var.items() :
            hoja.write(fila, 17, k, izquierda)
            hoja.write(fila, 18, v, derecha)
            fila += 1
        hoja.write(fila, 17, k, bajoIzq)
        hoja.write(fila, 18, v, bajoDrch)


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
    full.freeze_panes(1,0)
    #Comprobar si existen los archivos de variantes, leerlos y montarlos en el excel
    for a in arxius :
        if os.path.isfile(a) :
            aux = a.split(".")[0] # Recoger el nombre que tendra la pestana
            full = wb.add_worksheet(aux.capitalize())
            if aux == "cand" :
                full.activate()
            full.freeze_panes(1,0)
            dc = convertirArchivo(a)
            filaActual = crearCabecera(full, wb)
            filaActual = escribirVariantes(full, wb, dc, filaActual)
            ayudaPredictores(full, wb, filaActual+2)

    full = wb.add_worksheet("QC_stats")
    escribirEstadisticas(full, wb)
    print("INFO: Creado archivo excel, con nombre {}, con el resumen de resultados".format("{}.xlsx".format(nom)))
    wb.close()

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
