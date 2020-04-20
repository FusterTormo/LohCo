#!/usr/bin/python
# -*- coding: utf-8 -*-

import xlsxwriter
import getWebInfo
import sys
import os
import re
import time
import mmap

"""Constantes"""
minReads = 25 #Numero minimo de reads alterados para consderar una varinate
minVAF = 10 #Frecuencia minima para considerar una variante como candidata
#TODO editar documentacion!!!!!!!!!!!!!!!!
"""
Orden de las columnas para el archivo de texto.
El array esta dividido en bloques: definicion de la variante, informacion relevante del vcf, identificadores web en bases de datos de enfermedades,
MAF recogidas de ANNOVAR, MAF recogidas de myvariant.info, informacion de predictores recogida por ANNOVAR y datos crudos del vcf
"""
orden_all = ['Gene.refGene', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'ExonicFunc.refGene', 'GeneDetail.refGene', 'AAChange.refGene', "population_max", "predictor_summary", "Strand_bias_score",

'GT', 'GQ', 'DP', 'RFD', 'ALD', 'VAF', 'MT',

'avsnp150', 'cosmic70', 'CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG',

'gnomad_exome_AF', 'gnomad_exome_AF_popmax', 'gnomad_exome_AF_male', 'gnomad_exome_AF_female', 'gnomad_exome_AF_raw', 'gnomad_exome_AF_afr', 'gnomad_exome_AF_sas', 'gnomad_exome_AF_amr',
'gnomad_exome_AF_eas', 'gnomad_exome_AF_nfe', 'gnomad_exome_AF_fin', 'gnomad_exome_AF_asj', 'gnomad_exome_AF_oth', 'gnomad_exome_non_topmed_AF_popmax', 'gnomad_exome_non_neuro_AF_popmax',
'gnomad_exome_non_cancer_AF_popmax', 'gnomad_exome_controls_AF_popmax',
'gnomad_genome_AF', 'gnomad_genome_AF_popmax', 'gnomad_genome_AF_male', 'gnomad_genome_AF_female', 'gnomad_genome_AF_raw', 'gnomad_genome_AF_afr', 'gnomad_genome_AF_sas', 'gnomad_genome_AF_amr',
'gnomad_genome_AF_eas', 'gnomad_genome_AF_nfe', 'gnomad_genome_AF_fin', 'gnomad_genome_AF_asj', 'gnomad_genome_AF_oth', 'gnomad_genome_non_topmed_AF_popmax', 'gnomad_genome_non_neuro_AF_popmax',
'gnomad_genome_non_cancer_AF_popmax', 'controls_AF_popmax',
'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE','ExAC_OTH', 'ExAC_SAS',
'1000g2015aug_all', '1000g2015aug_afr', '1000g2015aug_amr', '1000g2015aug_eas', '1000g2015aug_eur', '1000g2015aug_sas',
'esp6500siv2_all', 'esp6500siv2_ea', 'esp6500siv2_aa',

'dbSNP_MAF',
'CADD_1000g_all', 'CADD_1000g_afr', 'CADD_1000g_amr', 'CADD_1000g_eur', 'CADD_1000g_eas', 'CADD_1000g_sas',
'dbNSFP_1000g_all', 'dbNSFP_1000g_afr', 'dbNSFP_1000g_amr', 'dbNSFP_1000g_eur', 'dbNSFP_1000g_eas', 'dbNSFP_1000g_sas',
'CADD_ESP6500_all', 'CADD_ESP6500_ea', 'CADD_ESP6500_aa',
'dbNSFP_esp6500_all', 'dbNSFP_esp6500_ea', 'dbNSFP_esp6500_aa',
'ExAC_ExAC_all', 'ExAC_ExAC_afr', 'ExAC_ExAC_amr', 'ExAC_ExAC_eas', 'ExAC_ExAC_fin', 'ExAC_ExAC_nfe', 'ExAC_ExAC_oth', 'ExAC_ExAC_sas',
'dbNSFP_ExAC_all', 'dbNSFP_ExAC_afr', 'dbNSFP_ExAC_amr', 'dbNSFP_ExAC_eas', 'dbNSFP_ExAC_fin', 'dbNSFP_ExAC_nfe', 'dbNSFP_ExAC_oth', 'dbNSFP_ExAC_sas',
'gNOMAD_Exome_all', 'gNOMAD_Exome_afr', 'gNOMAD_Exome_amr', 'gNOMAD_Exome_asj', 'gNOMAD_Exome_eas', 'gNOMAD_Exome_fin', 'gNOMAD_Exome_nfe', 'gNOMAD_Exome_oth', 'gNOMAD_Exome_popmax',
'gNOMAD_Exome_raw', 'gNOMAD_Exome_sas',
'gNOMAD_Genome_all', 'gNOMAD_Genome_afr', 'gNOMAD_Genome_amr', 'gNOMAD_Genome_asj', 'gNOMAD_Genome_eas', 'gNOMAD_Genome_fin', 'gNOMAD_Genome_nfe', 'gNOMAD_Genome_oth', 'gNOMAD_Genome_popmax',
'gNOMAD_Genome_raw',

'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score',
'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred',
'MutationAssessor_score', 'MutationAssessor_score_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore',
'PROVEAN_pred', 'VEST3_score', 'VEST3_rankscore', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'M-CAP_score', 'M-CAP_rankscore',
'M-CAP_pred', 'REVEL_score', 'REVEL_rankscore', 'MutPred_score', 'MutPred_rankscore', 'CADD_raw', 'CADD_raw_rankscore', 'CADD_phred', 'DANN_score', 'DANN_rankscore', 'fathmm-MKL_coding_score',
'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_pred', 'Eigen_coding_or_noncoding', 'Eigen-raw', 'Eigen-PC-raw', 'GenoCanyon_score', 'GenoCanyon_score_rankscore', 'integrated_fitCons_score',
'integrated_fitCons_score_rankscore', 'integrated_confidence_value', 'GERP++_RS', 'GERP++_RS_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP20way_mammalian',
'phyloP20way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore', 'SiPhy_29way_logOdds',
'SiPhy_29way_logOdds_rankscore', 'Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue',

'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE', 'TUMOR']


##############################
#    Orden definitivo        #
##############################

orden = ['muestra', 'IGV', 'IGV_analisis', 'Gene.refGene', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'ExonicFunc.refGene', 'GeneDetail.refGene', 'AAChange.refGene', "population_max",
"gNOMAD_Exome_popmax", "gNOMAD_Genome_popmax", "predictor_summary", "Strand_bias_score", 'GT', 'GQ', 'DP', 'RFD', 'ALD', 'VAF', 'MT', 'avsnp150', 'cosmic70', 'CLNALLELEID', 'CLNDN', 'CLNDISDB',
'CLNREVSTAT', 'CLNSIG']

#Columnas con informacion sobre MAF
mafCols = ['1000g2015aug_all', '1000g2015aug_afr', '1000g2015aug_amr', '1000g2015aug_eas', '1000g2015aug_eur', '1000g2015aug_sas', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN',
'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', 'gnomad_exome_AF', 'gnomad_exome_AF_popmax', 'gnomad_exome_AF_male', 'gnomad_exome_AF_female', 'gnomad_exome_AF_raw', 'gnomad_exome_AF_afr',
'gnomad_exome_AF_sas', 'gnomad_exome_AF_amr', 'gnomad_exome_AF_eas', 'gnomad_exome_AF_nfe', 'gnomad_exome_AF_fin', 'gnomad_exome_AF_asj', 'gnomad_exome_AF_oth',
'gnomad_exome_non_topmed_AF_popmax', 'gnomad_exome_non_neuro_AF_popmax','gnomad_exome_non_cancer_AF_popmax', 'gnomad_exome_controls_AF_popmax',
'gnomad_genome_AF', 'gnomad_genome_AF_popmax', 'gnomad_genome_AF_male', 'gnomad_genome_AF_female', 'gnomad_genome_AF_raw', 'gnomad_genome_AF_afr', 'gnomad_genome_AF_sas', 'gnomad_genome_AF_amr',
'gnomad_genome_AF_eas', 'gnomad_genome_AF_nfe', 'gnomad_genome_AF_fin', 'gnomad_genome_AF_asj', 'gnomad_genome_AF_oth', 'gnomad_genome_non_topmed_AF_popmax', 'gnomad_genome_non_neuro_AF_popmax',
'gnomad_genome_non_cancer_AF_popmax', 'controls_AF_popmax', 'esp6500siv2_all', 'esp6500siv2_ea', 'esp6500siv2_aa', 'dbSNP_MAF','CADD_1000g_all', 'CADD_1000g_afr', 'CADD_1000g_amr',
'CADD_1000g_eur', 'CADD_1000g_eas', 'CADD_1000g_sas', 'dbNSFP_1000g_all', 'dbNSFP_1000g_afr', 'dbNSFP_1000g_amr', 'dbNSFP_1000g_eur', 'dbNSFP_1000g_eas', 'dbNSFP_1000g_sas',
'CADD_ESP6500_all', 'CADD_ESP6500_ea', 'CADD_ESP6500_aa', 'dbNSFP_esp6500_all', 'dbNSFP_esp6500_ea', 'dbNSFP_esp6500_aa', 'ExAC_ExAC_all', 'ExAC_ExAC_afr', 'ExAC_ExAC_amr', 'ExAC_ExAC_eas',
'ExAC_ExAC_fin', 'ExAC_ExAC_nfe', 'ExAC_ExAC_oth', 'ExAC_ExAC_sas', 'dbNSFP_ExAC_all', 'dbNSFP_ExAC_afr', 'dbNSFP_ExAC_amr', 'dbNSFP_ExAC_eas', 'dbNSFP_ExAC_fin', 'dbNSFP_ExAC_nfe',
'dbNSFP_ExAC_oth', 'dbNSFP_ExAC_sas', 'gNOMAD_Exome_all', 'gNOMAD_Exome_afr', 'gNOMAD_Exome_amr', 'gNOMAD_Exome_asj', 'gNOMAD_Exome_eas', 'gNOMAD_Exome_fin', 'gNOMAD_Exome_nfe',
'gNOMAD_Exome_oth', 'gNOMAD_Exome_popmax', 'gNOMAD_Exome_raw', 'gNOMAD_Exome_sas', 'gNOMAD_Genome_all', 'gNOMAD_Genome_afr', 'gNOMAD_Genome_amr', 'gNOMAD_Genome_asj', 'gNOMAD_Genome_eas',
'gNOMAD_Genome_fin', 'gNOMAD_Genome_nfe', 'gNOMAD_Genome_oth', 'gNOMAD_Genome_popmax', 'gNOMAD_Genome_raw']

# Nombres de los genes guardados en el panel SMD actual
gens_panell = ["CSF3R", "MPL", "DNMT3A", "IDH1", "GATA2", "STAG1", "TET2", "NPM1", "EZH2", "RAD21", "JAK2", "SMC3", "ETV6", "PTPN11", "IDH2", "NF1", "CEBPA", "RUNX1", "ZRSR2",
    "BCOR", "GATA1", "STAG2", "BCORL1"]

def addFORMAT(dic, somatic) :
    """
    Separar los datos de las columnas INFO y FORMAT del vcf en columnas. Se devuelven los datos separados en un diccionario. Se comprueba si alguna clave esta ya presente en el diccionario
    """
    claves = dic["FORMAT"].split(":")
    valores = dic["SAMPLE"].split(":")
    if somatic == "somatico" :
        valores2 = dic["TUMOR"].split(":")
    minidic = {}
    i = 0
    #Separar los datos de la columna FORMAT en distintas columnas
    for c in claves :
        if c in dic.keys() :
            print "WARNING: Clau {} repetida en {}".format(c, dic.keys())
            sys.exit()
        else :
            minidic[c] = valores[i]
        i += 1
    #Separar los datos de la muestra somatic
    if somatic == "somatico" :
        i = 0
        for c in claves :
            newc = "{}_TUMOR".format(c)
            if newc in dic.keys() :
                print "WARNING: Clau {} repetida".format(newc)
            else :
                minidic[newc] = valores2[i]
            i += 1

    #Separar los datos de la columna INFO en distintas columnas
    datos = dic["INFO"].split(";")
    for d in datos :
        aux = d.split("=")
        aux[0] = "{}_INFO".format(aux[0])
        if aux[0] in dic.keys() or aux[0] in minidic.keys():
            print "WARNING: Clau {} repetida al muntar INFO"

        if len(aux) == 1 :
            minidic[aux[0]] = aux[0]
        else :
            minidic[aux[0]] = aux[1]

    return minidic

def storeData(path, som) :
    """Abrir un archivo con formato table_annovar y guardar todo el contenido en un diccionario"""
    cabecera = True
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

    print "INFO: Recollint informacio d'ANNOVAR"

    if som == "somatico" :
        claves.append("TUMOR")

    temp = {}
    dc = []
    stats = {} #Contar las variantes que tiene cada filtro
    with open(path, "r") as fi :
        for l in fi :
            if cabecera :
                cabecera = False
            else :
                aux = l.strip("\n").split("\t")
                if aux[144] in stats.keys() :
                    stats[aux[144]] += 1
                else :
                    stats[aux[144]] = 1

                if aux[144] == "PASS" : #NOTE solo las variantes que han pasado los filtros se guardan en el excel
                    i = 0
                    for c in claves :
                        temp[c] = aux[i]
                        i += 1

                    aux2 = addFORMAT(temp, som)
                    temp.update(aux2)
                    dc.append(temp)
                    temp = {}
    print "INFO: Resum dels filtres per la mostra"
    printStats(stats)
    print "INFO: Nomes es guardaran en l'excel aquelles variants que hagin passat els filtres de Strelka2"
    return dc

def addWebInfo(dc) :
    """Agregar la informacion de myvariant.info a aquellas variantes que sean exonicas"""
    stats = {}
    print "INFO: Recollint informacio des de myvariant.info"
    for l in dc :
        #Sacar estadisticas de los tipos de variantes que han pasado los filtros
        if l["Func.refGene"] in stats.keys() :
            stats[l["Func.refGene"]] += 1
        else :
            stats[l["Func.refGene"]] = 1
        if l["Func.refGene"] == "exonic" and l["ExonicFunc.refGene"] != "synonymous SNV" :
            aux = getWebInfo.getMyVariant(l["Chr"], l["Start"], l["Ref"], l["Alt"])
            l.update(aux)

    print "\nINFO: Resum de variants segons la seva posicio genomica"
    printStats(stats)

def write2Excel(dc, somatic, mostra="noName") :
    """Agregar al diccionario las columnas que se calculan y escribir toda la informacion ordenada en un archivo excel separando las variantes por pestanas segun pasen o no los filtros
    correspondientes del GESMD"""
    print "INFO: Filtrant i escrivint Excel"
    for l in dc :
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

    textName = "{}_excel.tsv".format(mostra)
    with open(textName, "w") as fi :
        fi.write("\t".join(orden_all))
        for f in dc :
            for k in orden_all :
                if k in f.keys() :
                    fi.write("{}\t".format(f[k]))
                else :
                    fi.write("NP\t")
            fi.write("\n")

    excelName = "{}.xlsx".format(mostra)
    tabs = 7
    completado = False
    while tabs > 0 and not completado :
        try :
            print "INFO: Guardant informacio en {} pestanyes".format(tabs)
            newExcel(excelName, mostra, somatic, dc, tabs)
            completado = True
        except UnicodeDecodeError :
            tabs -= 1
            print "WARNING: No es poden guardar totes les dades en la taula. Guardant {} pestanyes".format(tabs)

def newExcel(name, mostra, somatic, dc, tabs = 6) :
    if len(mostra) > 15 :
        mostra = mostra.split("_")[0]
    wb = xlsxwriter.Workbook(name, {"strings_to_numbers" : True, "tmpdir" : "/imppc/labs/solelab/ffuster2/tmp"})
    if tabs >= 1 :
        full = wb.add_worksheet("{}_Filtered_def".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)

    if tabs >= 2 :
        full = wb.add_worksheet("{}_CAND_panell".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        fila = 2
        for f in dc :
            if f["Gene.refGene"] in gens_panell and ((f["Func.refGene"] == "exonic" and f["ExonicFunc.refGene"] != "synonymous SNV") or f["Func.refGene"] == "splicing") :
                if f["population_max"] == "NA" or f["population_max"] < 0.01 :
                    if float(f["VAF"]) >= minVAF :
                        printLine(fila, f, full, wb)
                        fila += 1
        printHelp(fila, full, wb, somatic)

    if tabs >= 3 :
        full = wb.add_worksheet("{}_CAND".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        full.activate()
        fila = 2
        for f in dc :
            if (f["Func.refGene"] == "exonic" and f["ExonicFunc.refGene"] != "synonymous SNV") or f["Func.refGene"] == "splicing" :
                if f["population_max"] == "NA" or f["population_max"] < 0.01 :
                    if float(f["VAF"]) >= minVAF :
                        printLine(fila, f, full, wb)
                        fila += 1
        printHelp(fila, full, wb, somatic)

    if tabs >= 4 :
        full = wb.add_worksheet("{}_VAF<10%".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        fila = 2
        for f in dc :
            if (f["Func.refGene"] == "exonic" and f["ExonicFunc.refGene"] != "synonymous SNV") or f["Func.refGene"] == "splicing" :
                if f["population_max"] == "NA" or f["population_max"] < 0.01 :
                    if f["VAF"] < minVAF:
                        printLine(fila, f, full, wb)
                        fila += 1

    if tabs >= 5 :
        full = wb.add_worksheet("{}_MAF_>_0.01".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        fila = 2
        for f in dc :
            if (f["Func.refGene"] == "exonic" and f["ExonicFunc.refGene"] != "synonymous SNV") or f["Func.refGene"] == "splicing" :
                if f["population_max"] != "NA" and f["population_max"] >= 0.01 :
                    printLine(fila, f, full, wb)
                    fila += 1

    if tabs >= 6 :
        full = wb.add_worksheet("{}_conseq".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        fila = 2
        for f in dc :
            if (f["Func.refGene"] == "exonic" and f["ExonicFunc.refGene"] != "synonymous SNV") or f["Func.refGene"] == "splicing" :
                printLine(fila, f, full, wb)
                fila += 1

    if tabs >= 7 :
        full = wb.add_worksheet("{}_all_variants".format(mostra))
        printHeader(full, wb, somatic)
        full.freeze_panes(2,0)
        fila = 2
        for f in dc :
            printLine(fila, f, full, wb)
            fila += 1

    wb.close()

def printHelp(fila, hoja, libro, tipoVarCal) :
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

    fila += 2
    #NOTE: En la version actual de strelka2excel no se guardan los datos de los predictores en el excel. Esta informacion no es necesaria
    """
    #Ayuda para interpretar la prediccion de algunos de los predictores
    titol = "HELP for in silico predictors"
    sift = "SIFT\nD -> Deleterious\nT -> Benign"
    polyphen = "Polyphen\nD -> Probably damaging\nP -> Possibly damaging\nB -> Benign"
    lrt = "LRT\nD -> Deletereous\nN -> Neutral\nU -> Unknown"
    mtaster = "Mutation Taster\nA -> Disease causing automatic\nD -> Disease causing\nN -> Polymophism\nP -> Polymophism automatic"
    massessor = "Mutation Assessor\nH -> High (Deleterious)\nM -> Medium (Deleterious)\nL -> Low (Tolerated)\nN -> Neutral (Tolerated)"
    provean = "PROVEAN\nD -> Deleterious\nT -> Benign"
    fathmm = "FATHMM\nD -> Deleterious\nN -> Benign"
    msvm = "MetaSVM\nD -> Deleterious\nT -> Benign"
    mlr = "MetaLR\nD -> Deleterious\nT -> Benign"
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
    """
    #Ayuda para interpretar la columna summary_predictors
    titol = "HELP for 'predictor summary' column"
    body = "The column summarizes the prediction of SIFT,\nPolyphen2 HDIV, Polyphen2 HVAR, LRT, MutationTaster,\nMutationAssessor, FATHMM, PROVEAN, MetaSVM, and MetaLR.\n\n"
    body += "It enumerates the number of (D)eletereous, (T)olerated,\n and (U)nknown prediction\n\n"
    body += "So 2D, 7T, 1U\nmeans\n2 deleterious, 7 tolerated, and 1 unknown predictions."
    hoja.merge_range(fila, 10, fila, 16, titol, titulo)
    hoja.merge_range(fila+1, 10, fila+9, 16, body, bajo)
    #Ayuda para interpretar las calidades de la variante, dependiendo de si la muestra analizada es somatica o germinal
    titol = "HELP for quality column, depending on variant calling analysis"
    if tipoVarCal == "somatico" :
        body = "Quality column represents the probability to be a somatic variant.\nIf the frequency of the somatic variant is significantly different in the tumor than in the normal"
    else :
        body = "Quality column represents the score for variant sites."

    hoja.merge_range(fila+12, 10, fila+12, 18, titol, titulo)
    hoja.merge_range(fila+13, 10, fila+14, 18, body, bajo)

    #Ayuda para los coeficientes de strand-bias
    titol = "HELP for Strand bias ratio"
    if tipoVarCal == "germinal" :
        body = "Strand bias is calculated with the formula:\nreads_reference_forward + reads_alterated_forward - read_reference_reverse - reads_alterated_reverse"
    else :
        body = "If the variant is a SNV, the ratio is calculated getting the data in SNVSB column from VCF.\nIf the variant is an indel, the ratio cannot be calculated, due to lack of information"

    hoja.merge_range(fila+17, 10, fila+17, 19, titol, titulo)
    hoja.merge_range(fila+18, 10, fila+19, 19, body, bajo)

def printLine(fila, dic, hoja, excel) :
    """Crea estilos para la hoja excel y dibuja la fila"""
    rojo = excel.add_format({"bg_color" : "#FF4D4D"})
    verde = excel.add_format({"bg_color" : "#43F906"})
    amarillo = excel.add_format({"bg_color" : "#FFFF00"})
    naranja = excel.add_format({"bg_color" : "#FF8000"})
    columna = 2
    predictors = ["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred"]
    for o in orden :
        if o in dic.keys() :
            if o in mafCols or o == "population_max":
                if dic[o] != "NA" and dic[o] != "." and float(dic[o]) >= 0.01 :
                    hoja.write(fila, columna, dic[o], rojo)
                else :
                    hoja.write(fila, columna, dic[o])
            elif o in predictors : #Colorear la prediccion de los predictores in silico
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
        else :
            hoja.write(fila, columna, "NP")
        columna += 1

def printHeader(hoja, excel, somatic) :
    """Crear unos estilos para la hoja excel y dibujar la cabecera de la hoja con los estilos creados"""
    titulo = excel.add_format({'bold' : True,
        'align' : 'center',
        'border' : 1,
        'bg_color' :  '#B3E6FF',
        'font_size' : 13
    })
    #Titulos de las cabeceras
    superCabecera = [["Analysis & Interpretation",5], ["Variant description",6], ["Change description",4], ["Classification summary", 5], ["Quality params",7], ["dbSNP",1], ["COSMIC", 1],
    ["Clinvar",5]]

    cabecera = ["Analysis", "Interpretation", "Sample", "IGV", "IGV Analysis", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Function", "Exonic consequence", "Transcript", "AA change & transcript",
    "Max MAF in population DB", "Max MAF gNOMAD exome", "Max MAF gNOMAD genome", "Predictor summary", "Strand bias ratio", "Genome type", "Genome quality", "Depth", "Reference depth", "Alternative depth",
    "VAF", "Mutation type", "dbSNP", "COSMIC", "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG"]

    #Escribir las cabeceras en cada una de las pestanas
    cols = 0
    for l in superCabecera :
        if l[1] == 1 :
            hoja.write(0, cols, "{}".format(l[0]), titulo)
        else :
            hoja.merge_range(0, cols, 0, cols + l[1] -1, "{}".format(l[0]), titulo)
        cols += l[1]

    cols = 0
    for n in cabecera :
        hoja.write(1, cols, n, titulo)
        cols += 1

def printStats(dic) :
    count = 0
    for k, v in dic.iteritems() :
        print "\t{} -> {}".format(v, k)
        count += int(v)
    print "\t\t{} en total".format(count)

def somaORgerm(path) :
    cabecera = True
    queEs = ""
    with open(path) as fi :
        for l in fi :
            if cabecera :
                cabecera = False
            else :
                s = re.search("FDP", l)
                if s == None:
                    queEs = "germinal"
                else :
                    queEs = "somatico"
                break
    print "INFO: Comprovant si la mostra es somatica o germinal. Detectat {}".format(queEs)
    return queEs

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

arx = sys.argv[1]
if os.path.isfile(arx) :
    sg = somaORgerm(arx)
    dc = storeData(arx, sg)
    addWebInfo(dc)
    if len(sys.argv) > 2 :
        write2Excel(dc, sg, sys.argv[2])
    else :
        write2Excel(dc, sg)
