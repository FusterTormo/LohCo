#!/usr/bin/python
# -*- coding: utf-8 -*-

import xlsxwriter
import getWebInfo
import sys
import os
import re
import time


'''Constantes'''
maf = 0.01
#Sufijos de las hojas de Excel
sufijos = ["_Filtered_def", "_CAND", "_MAF_>_0.01", "_conseq", "_PASS_variants", "_all_variants"]
superCabecera = [["Analysis & Interpretation",5], ["Variant description",6], ["Quality params",7], ["Change description",4], ["dbSNP",1], ["ANNOVAR 1000g",6], ["ANNOVAR ExAC",8], ["ANNOVAR EVS",3],
["MAF dbSNP",1], ["CADD 1000g",6], ["dbNSFP 1000g",6], ["ExAC",8], ["dbNSFP ExAC",8], ["CADD EVS",3], ["dbNSFP EVS",3], ["Clinvar",5], ["COSMIC",1], ["Predictors",34], ["VCF additional info",4]]
cabecera = ["Analysis", "Interpretation", "Sample", "IGV", "IGV Analysis", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Mutation type", "Genome type", "Quality",
"Depth (Normal sample, Tumor sample if applicable)", "Reference depth", "Alt depth", "VAF", "Function", "Transcript", "Exonic consequence", "AA change & transcript", "avsnp147", "1000g2015aug_all",
"1000g2015aug_afr", "1000g2015aug_amr", "1000g2015aug_eas", "1000g2015aug_eur", "1000g2015aug_sas", "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
"esp6500siv2_all", "esp6500siv2_aa", "esp6500siv2_ea", "dbSNP_MAF", "CADD_1000g_all", "CADD_1000g_afr", "CADD_1000g_amr", "CADD_1000g_eas", "CADD_1000g_eur", "CADD_1000g_sas", "dbNSFP_1000g_all",
"dbNSFP_1000g_afr", "dbNSFP_1000g_amr",	"dbNSFP_1000g_eas", "dbNSFP_1000g_eur", "dbNSFP_1000g_sas", "ExAC_ExAC_all", "ExAC_ExAC_afr", "ExAC_ExAC_amr", "ExAC_ExAC_eas", "ExAC_ExAC_fin", "ExAC_ExAC_nfe",
"ExAC_ExAC_oth", "ExAC_ExAC_sas", "dbNSFP_ExAC_all", "dbNSFP_ExAC_afr", "dbNSFP_ExAC_amr", "dbNSFP_ExAC_eas", "dbNSFP_ExAC_fin", "dbNSFP_ExAC_nfe", "dbNSFP_ExAC_oth", "dbNSFP_ExAC_sas", "CADD_ESP6500_all",
"CADD_ESP6500_aa", "CADD_ESP6500_ea", "dbNSFP_esp6500_all", "dbNSFP_esp6500_aa", "dbNSFP_esp6500_ea", "CLINSIG", "CLNACC", "CLNDBN", "CLNDSDB", "CLNDSDBID", "cosmic70", "SIFT_score", "SIFT_pred",
"Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "LRT_score", "LRT_pred", "MutationTaster_score", "MutationTaster_pred", "MutationAssessor_score",
"MutationAssessor_pred", "FATHMM_score", "FATHMM_pred", "PROVEAN_score", "PROVEAN_pred", "VEST3_score", "CADD_raw", "CADD_phred", "DANN_score", "fathmm-MKL_coding_score", "fathmm-MKL_coding_pred",
"MetaSVM_score", "MetaSVM_pred", "MetaLR_score", "MetaLR_pred", "integrated_fitCons_score", "integrated_confidence_value", "GERP++_RS", "phyloP7way_vertebrate", "phyloP20way_mammalian",
"phastCons7way_vertebrate", "phastCons20way_mammalian", "SiPhy_29way_logOdds", "VCF FILTER", "INFO from VCF", "FORMAT from VCF", "FORMAT VALUES from VCF"]

def getTemps() :
    '''Formata l'hora actual del sistema i la retorna en format string
    @return String Hora actual en format HH:MM:SS'''
    t = time.localtime(time.time())
    return "{}:{}:{}".format(t.tm_hour, t.tm_min, t.tm_sec)

def extraureVAF(l, varCaller) :
  '''Devolver la frecuencia alelica de la variante que se pasa por parametro. Dependiendo del variant caller que se pase por parametro, los datos estaran situados en una u otra columna de la variante
  @param l {Tuple} variante de la cual se quiere extraer la frecuencia alelica de la variante
  @param varCaller {String} Tipo de variant caller en el que esta formateada la variante. Valores aceptados: strelka, lancet
  @return {List} (Coverage de la referencia, coverage del alterado, Frecuencia alelica de la variante en %)'''

  #Strelka sigue una formula para calcular el coverage que consiste, basicamente en ver el coverage que tiene cada una de las bases que corresponde a la variante y calcular la VAF usando dichos valoreS
  if varCaller == "strelka" :
      indexRef = -1
      indexAlt = -1
      index = -1
      formato = l[76].split(":")
      #Calcular coverage para SNV
      if len(l[3]) == 1 and len(l[4]) == 1 and l[3] != "-" and l[4] != "-" :
          ref = "{}U".format(l[3])
          alt = "{}U".format(l[4])
      else :
          ref = "TAR"
          alt = "TIR"

      for a in formato :
          index += 1
          if a == ref :
              indexRef = index
          if a == alt :
              indexAlt = index

      ref = l[78].split(":")[indexRef]
      alt = l[78].split(":")[indexAlt]
      covRef = ref.split(",")[0].strip()
      covAlt = alt.split(",")[0].strip()
  elif varCaller == "lancet" :
      aux = l[77].split(":")[1]
      covRef = aux.split(",")[0]
      covAlt = aux.split(",")[1]

  try :
      vaf = float(covAlt)/(float(covRef)+float(covAlt)) * 100
  except ZeroDivisionError :
      print "\tERROR: No es pot extraure la VAF amb les dades de  covAlt={}, covRef={}".format(covAlt, covRef)
      print "\n{}\n".format(l)
      vaf = 0

  return [covRef, covAlt, vaf]

def extraureCoverage(l, variantCaller) :
    '''Devolver el coverage que tiene la variante que se pasa por parametro. Para ello busca la columna DP dentro de la columna FORMAT. Dependiendo del variant caller, se extraera solo un valor (en caso
    de lancet), o dos, correspondientes a muestra normal y muestra tumor respectivamente (strelka)
    @param l {List} Variante de la cual se quiere extraer el coverage
    @param variantCaller Tipo de variant caller en el que esta formateada la variante. Valores aceptados: strelka, lancet
    @return {String} coverage de la variante'''
    formato = l[76].split(":")
    index = -1
    for a in formato :
        index += 1
        if a == "DP" :
            break

    if variantCaller == "strelka" :
        return "{}, {}".format(l[77].split(":")[index], l[78].split(":")[index])
    elif variantCaller == "lancet" :
        return "{}".format(l[77].split(":")[index])

def extraureQual(l, varCaller) :
    ''''Devolver la calidad genomica de la variante que se pasa por parametro. Dependiendo del variant caller que se pase, los datos estaran situados en una u otra columna
    @param l {List} variante de la cual se quiere extraer la calidad genomica
    @param varCaller {String} Tipo de variant caller en el que esta formateada la variante. Valores aceptados: strelka, lancet
    @return {Integer} Calidad genomica de la variante'''
    q = None
    if varCaller == "strelka" :
        #Parchazo por no reconocer la epxresion regunar "();QSS=\d+)"
        aux = l[75].split(";")[1]
        q = aux.split("=")[1]
    elif varCaller == "lancet" :
        q = l[73]

    return q

def extraureGenotip(l, variantCaller) :
    ''''Devolver el genotipo ("0/0","0/1", "1/1") de la variante que se pasa por parametro. Dependiendo del variant caller que se pase, el dato estara situado en una u otra columna
    @param l {List} variante de la cual se quiere extraer el genotipo
    @param varCaller {String} Tipo de variant caller en el que esta formateada la variante. Valores aceptados: strelka, lancet
    @return {String} Genotipo de la variante. Posibles valores "0/0" (referencia), "0/1" (heterocigota), "1/1" (homocigota)'''
    geno = None

    if variantCaller == "strelka" :
        geno = re.search('(NT=(\w){3})', l[75])
        if geno.group(0) == "NT=ref" :
            geno = "0/0"
        elif geno.group(0) == "NT=het" :
            geno = "0/1"
        elif geno.group(0) == "NT=hom" :
             geno = "1/1"
    elif variantCaller == "lancet" :
        geno = l[77].split(":")[0]

    return geno

def extraureFlag(l, variantCaller) :
    '''Devolver el parametro de donde se ha encontrado la variante. Dependiendo del variant caller, este dato estara en una columna u otra. Se puede devolver:
            SOMATIC en caso de que la variante solo este en la muestra tumoral
            SHARED en caso de que la variante este en ambas muestras (control y tumoral)º
            NORMAL en caso de que la variante este en la muestra control y no este en la muestra tumoral
    @param l {Tuple} variante de la cual se quiere extraer el parametro
    @param varCaller {String} Tipo de variant caller en el que esta formateada la variante. Valores aceptados: strelka, lancet
    @return {String} Genotipo de la variante. Posibles valores "SOMATIC", "SHARED", "NORMAL"'''
    flag = ""
    if variantCaller == "strelka" :
        if l[75].find("SOMATIC") > -1 :
            flag = "SOMATIC"
        elif l[75].find("SHARED") > -1 :
            flag = "SHARED"
        elif l[75].find("NORMAL") > -1 :
            flag = "NORMAL"
        else :
            flag = "NA"
    elif variantCaller == "lancet" :
        if l[75].find("SOMATIC") > -1 :
            flag = "SOMATIC"
        elif l[75].find("SHARED") > -1 :
            flag = "SHARED"
        elif l[75].find("NORMAL") > -1 :
            flag = "NORMAL"
        else :
            flag = "NA"

    return flag

def extraureFiltre(l) :
    return l[74]

def combinaOrdena(l, w, varCaller, i) :
    '''Juntar los datos locales y web de la variante. Adaptar los datos  al formato tipico que pongo en los excel.
    @param {List} Datos de la variante extraidos del fichero
    @parm {Dict} Datos de la variante extraidos de myvariant.info. En caso de no encontrarse informacion, este parametro tendra valor None
    @param {String} identificador del paciente. Sera el primer dato agregue a la lista de retorno
    @return {List} Variante en el formato especifico para el excel'''
    #Columnas de Variant description: Chr, ini, fin, ref, alt
    rt = [i, l[6], l[0], l[1], l[2], l[3], l[4]]

    #Columnas de Quality params: tipo de variante (SOMATIC, SHARED o NORMAL), genome_type, genome_quality, total_depth, reference_depth, alt_depth y VAF
    rt.append(extraureFlag(l,varCaller))
    rt.append(extraureGenotip(l,varCaller))
    rt.append(extraureQual(l,varCaller))
    rt.append(extraureCoverage(l, varCaller))
    rt = rt + extraureVAF(l, varCaller)

    #Columnas de change description
    rt.append(l[5])
    rt.append(l[7])
    rt.append(l[8])
    rt.append(l[9])

    #Parche para evitar puntos en la anotacion de ANNOVAR
    for it in range(10, 28) :
        if l[it] == "." :
            l[it] = "NA"

    #Columnas de bases de datos poblacionales locales
    rt.append(l[10]) #dbSNP147
    rt.append(l[11]) #1000genomes
    rt.append(l[12])
    rt.append(l[13])
    rt.append(l[14])
    rt.append(l[15])
    rt.append(l[16])
    rt.append(l[17]) #ExAC
    rt.append(l[18])
    rt.append(l[19])
    rt.append(l[20])
    rt.append(l[12])
    rt.append(l[22])
    rt.append(l[23])
    rt.append(l[24])
    rt.append(l[25]) #Exome Variant Server
    rt.append(l[26])
    rt.append(l[27])

    if w != None :
        #Columnas de bases de datos poblacionales de la web (myvariant.info)
        rt.append(w['dbSNP_MAF']) #dbSNP
        rt.append(w['CADD_1000g_all']) #CADD 1000genomes
        rt.append(w['CADD_1000g_afr'])
        rt.append(w['CADD_1000g_amr'])
        rt.append(w['CADD_1000g_eas'])
        rt.append(w['CADD_1000g_eur'])
        rt.append(w['CADD_1000g_sas'])
        rt.append(w['dbNSFP_1000g_all']) #dbnsfp 1000genomes
        rt.append(w['dbNSFP_1000g_afr'])
        rt.append(w['dbNSFP_1000g_amr'])
        rt.append(w['dbNSFP_1000g_eas'])
        rt.append(w['dbNSFP_1000g_eur'])
        rt.append(w['dbNSFP_1000g_sas'])
        rt.append(w['ExAC_ExAC_all']) #ExAC de la web
        rt.append(w['ExAC_ExAC_afr'])
        rt.append(w['ExAC_ExAC_amr'])
        rt.append(w['ExAC_ExAC_eas'])
        rt.append(w['ExAC_ExAC_fin'])
        rt.append(w['ExAC_ExAC_nfe'])
        rt.append(w['ExAC_ExAC_oth'])
        rt.append(w['ExAC_ExAC_sas'])
        rt.append(w['dbNSFP_ExAC_all']) # #dbNSFP ExAC
        rt.append(w['dbNSFP_ExAC_afr'])
        rt.append(w['dbNSFP_ExAC_amr'])
        rt.append(w['dbNSFP_ExAC_eas'])
        rt.append(w['dbNSFP_ExAC_fin'])
        rt.append(w['dbNSFP_ExAC_nfe'])
        rt.append(w['dbNSFP_ExAC_oth'])
        rt.append(w['dbNSFP_ExAC_sas'])
        rt.append(w['CADD_ESP6500_all']) #CADD EVS
        rt.append(w['CADD_ESP6500_aa'])
        rt.append(w['CADD_ESP6500_ea'])
        rt.append(w['dbNSFP_esp6500_all']) #dbNSFP EVS
        rt.append(w['dbNSFP_esp6500_aa'])
        rt.append(w['dbNSFP_esp6500_ea'])
    else :
        na = 'NA'
        rt = rt + [na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na, na]

    #Columnas de COSMIC y ClinVar
    rt.append(l[29])
    rt.append(l[30])
    rt.append(l[31])
    rt.append(l[32])
    rt.append(l[33])

    #Columnas restantes que contienen los predictores
    rt.append(l[33])
    rt.append(l[34])
    rt.append(l[35])
    rt.append(l[36])
    rt.append(l[37])
    rt.append(l[38])
    rt.append(l[39])
    rt.append(l[40])
    rt.append(l[41])
    rt.append(l[42])
    rt.append(l[43])
    rt.append(l[44])
    rt.append(l[45])
    rt.append(l[46])
    rt.append(l[47])
    rt.append(l[48])
    rt.append(l[49])
    rt.append(l[50])
    rt.append(l[51])
    rt.append(l[52])
    rt.append(l[53])
    rt.append(l[54])
    rt.append(l[55])
    rt.append(l[56])
    rt.append(l[57])
    rt.append(l[58])
    rt.append(l[59])
    rt.append(l[60])
    rt.append(l[61])
    rt.append(l[62])
    rt.append(l[63])
    rt.append(l[64])
    rt.append(l[65])
    rt.append(l[66])
    rt.append(l[67])

    #Columnas recogidas del vcf. Solo pongo las columnas correspondientes la columna filtro, columna info y las columas format
    rt.append(l[74])
    rt.append(l[75])
    rt.append(l[76])
    rt.append(l[77])
    if len(l) == 79 :
        rt.append(l[78])

    return rt

def escribirExcel(x,i,excel) :
    '''Escribir la variante en la pestana Excel correspondiente de acuerdo a las especificaciones de cada una de las pestanas
    Pestanas (orden de sufijos):
            5. Todas las variantes reportadas
            4. Variantes que han pasado todos los filtros del variant caller (columna FILTER = PASS)
            3. Variantes que cumplen 2 y son splicing o exonicas (no sinonimas)
            2. Variantes que cumplen 3 pero tienen alguna de las MAF >= 1%
            1. Variantes que cumplen 3 y todas las MAF son < 1% o NA
            0. Hoja vacia donde guardar la lista de candidatas finales'''

    '''Crear sub listas de la lista principal con los datos de cada una de las variantes.
    Escribir cada una de las listas en su hoja correspondiente de Excel.
    Crear la hoja en blanco para las variantes candidatas'''
    subl1 = []
    subl2 = []
    subl3 = []
    subl4 = []
    print "INFO\t{}\tSeparant les dades segons filtres".format(getTemps())
    for l in excel :
        if l[111] == "PASS" :
            subl4.append(l)
            if l[14] == "exonic" and l[16] != "synonymous SNV" :
                subl3.append(l)

                filtroMaf = True
                for it in range(19,71) :
                    if l[it] != 'NA' and float(l[it]) >= maf :
                        filtroMaf = False
                        break

                if filtroMaf :
                    subl1.append(l)
                else :
                    subl2.append(l)

    print "INFO\t{}\tGuardant les dades en Excel".format(getTemps())

    wb = xlsxwriter.Workbook("{}.xlsx".format(x), {"strings_to_numbers" : True})
    full = wb.add_worksheet("{}{}".format(i, sufijos[0])) #Hoja vacio
    crearCabecera(full, wb)
    full.freeze_panes(2,0)

    full = wb.add_worksheet("{}{}".format(i, sufijos[1])) #Variantes candidatas
    full.activate()
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    pintarLinea(subl1, full, wb)

    full = wb.add_worksheet("{}{}".format(i, sufijos[2]))
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    pintarLinea(subl2, full, wb)

    full = wb.add_worksheet("{}{}".format(i, sufijos[3]))
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    pintarLinea(subl3, full, wb)

    full = wb.add_worksheet("{}{}".format(i, sufijos[4]))
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    pintarLinea(subl4, full, wb)

    full = wb.add_worksheet("{}{}".format(i, sufijos[5]))
    crearCabecera(full, wb)
    full.freeze_panes(2,0)
    pintarLinea(excel, full, wb)

    wb.close()

def escribirTsv(i, lista) :
    '''Escribir las variantes en los ficheros separados por tablas correspondientes'''

    for j in range(1,6) :
        with open("{}{}.tsv".format(i, sufijos[j]), "w") as fi :
            for tit in superCabecera :
                fi.write(tit[0])
                fi.write("\t"*tit[1])
                fi.write("\n")

            fi.write("\t".join(cabecera))
            fi.write("\n")

    for l in lista :
        chars = [str(x) for x in l]
        with open("{}{}.tsv".format(i,sufijos[5]), "a") as fi:
            fi.write("\t".join(chars))
            fi.write("\n")

        if l[111] == "PASS" :
            with open("{}{}.tsv".format(i, sufijos[4]), "a") as fi :
                fi.write("\t".join(chars))
                fi.write("\n")

            if l[14] == "exonic" and l[16] != "synonymous SNV" :
                with open("{}{}.tsv".format(i, sufijos[3]), "a") as fi :
                    fi.write("\t".join(chars))
                    fi.write("\n")

                filtroMaf = True
                for it in range(19, 71) :
                    if l[it] != 'NA' and float(l[it]) >= maf :
                        filtroMaf = False
                        break

                if filtroMaf :
                    with open("{}{}.tsv".format(i, sufijos[1]), "a") as fi :
                        fi.write("\t".join(chars))
                        fi.write("\n")
                else :
                    with open("{}{}.tsv".format(i, sufijos[2]), "a") as fi :
                        fi.write("\t".join(chars))
                        fi.write("\n")

def pintarLinea(lista,ws, wb) :
    '''Escribir las variantes dentro de la hoja pasada por parametro'''
    #Estilos para marcar las celdas
    fondo_rojo = wb.add_format({ 'bg_color' : '#FF9999'})
    fondo_naranja = wb.add_format({ 'bg_color' : '#FFD9B3'})
    fondo_verde = wb.add_format({ 'bg_color' : '#99FF99'})
    fondo_amarillo = wb.add_format({ 'bg_color' : '#FFFF00'})
    fondo_azul = wb.add_format({'bold' : True,
        'align' : 'center',
        'border' : 1,
        'bg_color' :  '#B3E6FF',
        'font_size' : 13
    })

    rows = 2
    for l in lista :
        ws.write(rows, 0, "") #Analysis
        ws.write(rows, 1, "") #Interpreattion
        ws.write(rows, 2, "{}".format(l[0])) #Sample
        if len(lista) < 1000 :
            enlace = "http://localhost:60151/goto?locus={}:{}".format(l[2], l[3])
        else :
            enlace = "localhost:60151/goto?locus={}:{}".format(l[2], l[3])
        ws.write(rows, 3, enlace) #IGV
        ws.write(rows, 4, "") #IGV Analysis
        ws.write(rows, 5, "{}".format(l[1])) #Gene
        ws.write(rows, 6, "{}".format(l[2])) #Chr
        ws.write(rows, 7, "{}".format(l[3])) #Start
        ws.write(rows, 8, "{}".format(l[4])) #End
        ws.write(rows, 9, "{}".format(l[5])) #Ref
        ws.write(rows, 10, "{}".format(l[6])) #Alt
        ws.write(rows, 11, "{}".format(l[7])) #Mutation type
        ws.write(rows, 12, "{}".format(l[8])) #Genome type
        ws.write(rows, 13, "{}".format(l[9])) #Quality
        ws.write(rows, 14, "{}".format(l[10])) #Depth
        ws.write(rows, 15, "{}".format(l[11])) #Ref depth
        ws.write(rows, 16, "{}".format(l[12])) #Alt depth
        ws.write(rows, 17, "{}".format(l[13])) #VAF
        ws.write(rows, 18, "{}".format(l[14])) #Function
        ws.write(rows, 19, "{}".format(l[15])) #Transcript
        ws.write(rows, 20, "{}".format(l[16])) #Exonic consequence
        ws.write(rows, 21, "{}".format(l[17])) #AA change
        ws.write(rows, 22, "{}".format(l[18])) #dbSNP id

        if l[19] == 'NA' or float(l[19]) < maf :
            ws.write(rows, 23, "{}".format(l[19])) #1000g2015aug_all
        else :
            ws.write(rows, 23, "{}".format(l[19]), fondo_rojo)

        if l[20] == 'NA' or float(l[20]) < maf :
            ws.write(rows, 24, "{}".format(l[20])) #1000g2015aug_afr
        else :
            ws.write(rows, 24, "{}".format(l[20]), fondo_rojo)
        if l[21] == 'NA' or float(l[21]) < maf :
            ws.write(rows, 25, "{}".format(l[21])) #1000g2015aug_amr
        else :
            ws.write(rows, 25, "{}".format(l[21]), fondo_rojo)
        if l[22] == 'NA' or float(l[22]) < maf :
            ws.write(rows, 26, "{}".format(l[22])) #1000g2015aug_eas
        else :
            ws.write(rows, 26, "{}".format(l[22]), fondo_rojo)
        if l[23] == 'NA' or float(l[23]) < maf :
            ws.write(rows, 27, "{}".format(l[23])) #1000g2015aug_eur
        else :
            ws.write(rows, 27, "{}".format(l[23]), fondo_rojo)
        if l[24] == 'NA' or float(l[24]) < maf :
            ws.write(rows, 28, "{}".format(l[24])) #1000g2015aug_sas
        else :
            ws.write(rows, 28, "{}".format(l[24]), fondo_rojo)
        if l[25] == 'NA' or float(l[25]) < maf :
            ws.write(rows, 29, "{}".format(l[25])) #ExAC_ALL
        else :
            ws.write(rows, 29, "{}".format(l[25]), fondo_rojo)
        if l[26] == 'NA' or float(l[26]) < maf :
            ws.write(rows, 30, "{}".format(l[26])) #ExAC_AFR
        else :
            ws.write(rows, 30, "{}".format(l[26]), fondo_rojo)
        if l[27] == 'NA' or float(l[27]) < maf :
            ws.write(rows, 31, "{}".format(l[27])) #ExAC_AMR
        else :
            ws.write(rows, 31, "{}".format(l[27]), fondo_rojo)
        if l[28] == 'NA' or float(l[28]) < maf :
            ws.write(rows, 32, "{}".format(l[28])) #ExAC_EAS
        else :
            ws.write(rows, 32, "{}".format(l[28]), fondo_rojo)
        if l[29] == 'NA' or float(l[29]) < maf :
            ws.write(rows, 33, "{}".format(l[29])) #ExAC_FIN
        else :
            ws.write(rows, 33, "{}".format(l[29]), fondo_rojo)
        if l[30] == 'NA' or float(l[30]) < maf :
            ws.write(rows, 34, "{}".format(l[30])) #ExAC_NFE
        else :
            ws.write(rows, 34, "{}".format(l[30]), fondo_rojo)
        if l[31] == 'NA' or float(l[31]) < maf :
            ws.write(rows, 35, "{}".format(l[31])) #ExAC_OTH
        else :
            ws.write(rows, 35, "{}".format(l[31]), fondo_rojo)
        if l[32] == 'NA' or float(l[32]) < maf :
            ws.write(rows, 36, "{}".format(l[32])) #ExAC_SAS
        else :
            ws.write(rows, 36, "{}".format(l[32]), fondo_rojo)
        if l[33] == 'NA' or float(l[33]) < maf :
            ws.write(rows, 37, "{}".format(l[33])) #esp6500siv2_all
        else :
            ws.write(rows, 37, "{}".format(l[33]), fondo_rojo)
        if l[34] == 'NA' or float(l[34]) < maf :
            ws.write(rows, 38, "{}".format(l[34])) #esp6500siv2_aa
        else :
            ws.write(rows, 38, "{}".format(l[34]), fondo_rojo)
        if l[35] == 'NA' or float(l[35]) < maf :
            ws.write(rows, 39, "{}".format(l[35])) #esp6500siv2_ea
        else :
            ws.write(rows, 39, "{}".format(l[35]), fondo_rojo)
        if l[36] == 'NA' or float(l[36]) < maf :
            ws.write(rows, 40, "{}".format(l[36])) #dbSNP_MAF
        else :
            ws.write(rows, 40, "{}".format(l[36]), fondo_rojo)
        if l[37] == 'NA' or float(l[37]) < maf :
            ws.write(rows, 41, "{}".format(l[37])) #CADD_1000g_all
        else :
            ws.write(rows, 41, "{}".format(l[37]), fondo_rojo)
        if l[38] == 'NA' or float(l[38]) < maf :
            ws.write(rows, 42, "{}".format(l[38])) #CADD_1000g_afr
        else :
            ws.write(rows, 42, "{}".format(l[38]), fondo_rojo)
        if l[39] == 'NA' or float(l[39]) < maf :
            ws.write(rows, 43, "{}".format(l[39])) ##CADD_1000g_amr
        else :
            ws.write(rows, 43, "{}".format(l[39]), fondo_rojo)
        if l[40] == 'NA' or float(l[40]) < maf :
            ws.write(rows, 44, "{}".format(l[40])) #CADD_1000g_eas
        else :
            ws.write(rows, 44, "{}".format(l[40]), fondo_rojo)
        if l[41] == 'NA' or float(l[41]) < maf :
            ws.write(rows, 45, "{}".format(l[41])) #CADD_1000g_eur
        else :
            ws.write(rows, 45, "{}".format(l[41]), fondo_rojo)
        if l[42] == 'NA' or float(l[42]) < maf :
            ws.write(rows, 46, "{}".format(l[42])) #CADD_1000g_sas
        else :
            ws.write(rows, 46, "{}".format(l[42]), fondo_rojo)
        if l[43] == 'NA' or float(l[43]) < maf :
            ws.write(rows, 47, "{}".format(l[43])) #dbNSFP_1000g_all
        else :
            ws.write(rows, 47, "{}".format(l[43]), fondo_rojo)
        if l[44] == 'NA' or float(l[44]) < maf :
            ws.write(rows, 48, "{}".format(l[44])) #dbNSFP_1000g_afr
        else :
            ws.write(rows, 48, "{}".format(l[44]), fondo_rojo)
        if l[45] == 'NA' or float(l[45]) < maf :
            ws.write(rows, 49, "{}".format(l[45])) #dbNSFP_1000g_amr
        else :
            ws.write(rows, 49, "{}".format(l[45]), fondo_rojo)
        if l[46] == 'NA' or float(l[46]) < maf :
            ws.write(rows, 50, "{}".format(l[46])) #dbNSFP_1000g_eas
        else :
            ws.write(rows, 50, "{}".format(l[46]), fondo_rojo)
        if l[47] == 'NA' or float(l[47]) < maf :
            ws.write(rows, 51, "{}".format(l[47])) #dbNSFP_1000g_eur
        else :
            ws.write(rows, 51, "{}".format(l[47]), fondo_rojo)
        if l[48] == 'NA' or float(l[48]) < maf :
            ws.write(rows, 52, "{}".format(l[48])) #dbNSFP_1000g_sas
        else :
            ws.write(rows, 52, "{}".format(l[48]), fondo_rojo)
        if l[49] == 'NA' or float(l[49]) < maf :
            ws.write(rows, 53, "{}".format(l[49])) #ExAC_ExAC_all
        else :
            ws.write(rows, 53, "{}".format(l[49]), fondo_rojo)
        if l[50] == 'NA' or float(l[50]) < maf :
            ws.write(rows, 54, "{}".format(l[50])) #ExAC_ExAC_afr
        else :
            ws.write(rows, 54, "{}".format(l[50]), fondo_rojo)
        if l[51] == 'NA' or float(l[51]) < maf :
            ws.write(rows, 55, "{}".format(l[51])) #ExAC_ExAC_amr
        else :
            ws.write(rows, 55, "{}".format(l[51]), fondo_rojo)
        if l[52] == 'NA' or float(l[52]) < maf :
            ws.write(rows, 56, "{}".format(l[52])) #ExAC_ExAC_eas
        else :
            ws.write(rows, 56, "{}".format(l[52]), fondo_rojo)
        if l[53] == 'NA' or float(l[53]) < maf :
            ws.write(rows, 57, "{}".format(l[53])) #ExAC_ExAC_fin
        else :
            ws.write(rows, 57, "{}".format(l[53]), fondo_rojo)
        if l[54] == 'NA' or float(l[54]) < maf :
            ws.write(rows, 58, "{}".format(l[54])) #ExAC_ExAC_nfe
        else :
            ws.write(rows, 58, "{}".format(l[54]), fondo_rojo)
        if l[55] == 'NA' or float(l[55]) < maf :
            ws.write(rows, 59, "{}".format(l[55])) #ExAC_ExAC_oth
        else :
            ws.write(rows, 59, "{}".format(l[55]), fondo_rojo)
        if l[56] == 'NA' or float(l[56]) < maf :
            ws.write(rows, 60, "{}".format(l[56])) #ExAC_ExAC_sas
        else :
            ws.write(rows, 60, "{}".format(l[56]), fondo_rojo)
        if l[57] == 'NA' or float(l[57]) < maf :
            ws.write(rows, 61, "{}".format(l[57])) #dbNSFP_ExAC_all
        else :
            ws.write(rows, 61, "{}".format(l[57]), fondo_rojo)
        if l[58] == 'NA' or float(l[58]) < maf :
            ws.write(rows, 62, "{}".format(l[58])) #dbNSFP_ExAC_afr
        else :
            ws.write(rows, 62, "{}".format(l[58]), fondo_rojo)
        if l[59] == 'NA' or float(l[59]) < maf :
            ws.write(rows, 63, "{}".format(l[59])) #dbNSFP_ExAC_amr
        else :
            ws.write(rows, 63, "{}".format(l[59]), fondo_rojo)
        if l[60] == 'NA' or float(l[60]) < maf :
            ws.write(rows, 64, "{}".format(l[60])) #dbNSFP_ExAC_eas
        else :
            ws.write(rows, 64, "{}".format(l[60]), fondo_rojo)
        if l[61] == 'NA' or float(l[61]) < maf :
            ws.write(rows, 65, "{}".format(l[61])) #dbNSFP_ExAC_fin"
        else :
            ws.write(rows, 65, "{}".format(l[61]), fondo_rojo)
        if l[62] == 'NA' or float(l[62]) < maf :
            ws.write(rows, 66, "{}".format(l[62])) #dbNSFP_ExAC_nfe
        else :
            ws.write(rows, 66, "{}".format(l[62]), fondo_rojo)
        if l[63] == 'NA' or float(l[63]) < maf :
            ws.write(rows, 67, "{}".format(l[63])) #dbNSFP_ExAC_oth
        else :
            ws.write(rows, 67, "{}".format(l[63]), fondo_rojo)
        if l[64] == 'NA' or float(l[64]) < maf :
            ws.write(rows, 68, "{}".format(l[64])) #dbNSFP_ExAC_sas
        else :
            ws.write(rows, 68, "{}".format(l[64]), fondo_rojo)
        if l[65] == 'NA' or float(l[65]) < maf :
            ws.write(rows, 69, "{}".format(l[65])) #CADD_ESP6500_all
        else :
            ws.write(rows, 69, "{}".format(l[65]), fondo_rojo)
        if l[66] == 'NA' or float(l[66]) < maf :
            ws.write(rows, 70, "{}".format(l[66])) #CADD_ESP6500_aa
        else :
            ws.write(rows, 70, "{}".format(l[66]), fondo_rojo)
        if l[67] == 'NA' or float(l[67]) < maf :
            ws.write(rows, 71, "{}".format(l[67])) #CADD_ESP6500_ea
        else :
            ws.write(rows, 71, "{}".format(l[67]), fondo_rojo)
        if l[68] == 'NA' or float(l[68]) < maf :
            ws.write(rows, 72, "{}".format(l[68])) #dbNSFP_esp6500_all
        else :
            ws.write(rows, 72, "{}".format(l[68]), fondo_rojo)
        if l[69] == 'NA' or float(l[69]) < maf :
            ws.write(rows, 73, "{}".format(l[69])) #dbNSFP_esp6500_aa
        else :
            ws.write(rows, 73, "{}".format(l[69]), fondo_rojo)
        if l[70] == 'NA' or float(l[70]) < maf :
            ws.write(rows, 74, "{}".format(l[70])) #dbNSFP_esp6500_ea
        else :
            ws.write(rows, 74, "{}".format(l[70]), fondo_rojo)

        #Columnas identificadores de Clinvar y COSMIC
        ws.write(rows, 75, "{}".format(l[71])) #CLINSIG
        ws.write(rows, 76, "{}".format(l[72])) #CLNACC
        ws.write(rows, 77, "{}".format(l[73])) #CLNDBN
        ws.write(rows, 78, "{}".format(l[74])) #CLNDSDB
        ws.write(rows, 79, "{}".format(l[75])) #CLNDSDBID
        ws.write(rows, 80, "{}".format(l[76]))#cosmic70

        #Predictores
        ws.write(rows, 81, "{}".format(l[77])) #SIFT_score
        if l[78] == "D" : #SIFT_pred
            ws.write(rows, 82, "{}".format(l[78]), fondo_rojo)
        elif l[78] == "T" :
            ws.write(rows, 82, "{}".format(l[78]), fondo_verde)
        else :
            ws.write(rows, 82, "{}".format(l[78]))
        ws.write(rows, 83, "{}".format(l[79])) #Polyphen2_HDIV_score
        if l[80] == "D" : #Polyphen2_HDIV_pred
            ws.write(rows, 84, "{}".format(l[80]), fondo_rojo)
        elif l[80] == "P" :
            ws.write(rows, 84, "{}".format(l[80]), fondo_naranja)
        elif l[80] == "B" :
            ws.write(rows, 84, "{}".format(l[80]), fondo_verde)
        else :
            ws.write(rows, 84, "{}".format(l[80]))
        ws.write(rows, 85, "{}".format(l[81])) #Polyphen2_HVAR_score
        if l[82] == "D" : #Polyphen2_HVAR_pred
            ws.write(rows, 86, "{}".format(l[82]), fondo_rojo)
        elif l[82] == "P" :
            ws.write(rows, 86, "{}".format(l[82]), fondo_naranja)
        elif l[82] == "B" :
            ws.write(rows, 86, "{}".format(l[82]), fondo_verde)
        else :
            ws.write(rows, 86, "{}".format(l[82]))
        ws.write(rows, 87, "{}".format(l[83])) #LRT_score
        if l[84] == "D" : #LRT_pred
            ws.write(rows, 88, "{}".format(l[84]), fondo_rojo)
        elif l[84] == "U" :
            ws.write(rows, 88, "{}".format(l[84]), fondo_amarillo)
        elif l[84] == "N":
            ws.write(rows, 88, "{}".format(l[84]), fondo_verde)
        else :
            ws.write(rows, 88, "{}".format(l[84]))
        ws.write(rows, 89, "{}".format(l[85])) #MutationTaster_score
        if l[86] == "A" or l[86] == "D" : #MutationTaster_pred
            ws.write(rows, 90, "{}".format(l[86]), fondo_rojo)
        elif l[86] == "N" or l[86] == "P" :
            ws.write(rows, 90, "{}".format(l[86]), fondo_verde)
        else :
            ws.write(rows, 90, "{}".format(l[86]))
        ws.write(rows, 91, "{}".format(l[87])) #MutationAssessor_score
        if l[88] == "H" : #MutationAssessor_pred
            ws.write(rows, 92, "{}".format(l[88]), fondo_rojo)
        elif l[88] == "M" :
            ws.write(rows, 92, "{}".format(l[88]), fondo_naranja)
        elif l[88] == "L" or l[88] == "N" :
            ws.write(rows, 92, "{}".format(l[88]), fondo_verde)
        else :
            ws.write(rows, 92, "{}".format(l[88]))
        ws.write(rows, 93, "{}".format(l[89])) #FATHMM_score
        if l[90] == "D" : #FATHMM_pred
            ws.write(rows, 94, "{}".format(l[90]), fondo_rojo)
        elif l[90] == "T" :
            ws.write(rows, 94, "{}".format(l[90]), fondo_verde)
        else :
            ws.write(rows, 94, "{}".format(l[90]))
        ws.write(rows, 95, "{}".format(l[91])) #PROVEAN_score
        if l[92] == "D" : #PROVEAN_pred
            ws.write(rows, 96, "{}".format(l[92]), fondo_rojo)
        elif l[92] == "T" :
            ws.write(rows, 96, "{}".format(l[92]), fondo_verde)
        else :
            ws.write(rows, 96, "{}".format(l[92]))
        ws.write(rows, 97, "{}".format(l[93])) #VEST3_score
        ws.write(rows, 98, "{}".format(l[94])) #CADD_raw
        ws.write(rows, 99, "{}".format(l[95])) #CADD_phred
        ws.write(rows, 100, "{}".format(l[96])) #DANN_score
        ws.write(rows, 101, "{}".format(l[97])) #fathmm-MKL_coding_score
        if l[98] == "D" : #fathmm-MKL_coding_pred
            ws.write(rows, 102, "{}".format(l[98]), fondo_rojo)
        elif l[98] == "N" :
            ws.write(rows, 102, "{}".format(l[98]), fondo_verde)
        else :
            ws.write(rows, 102, "{}".format(l[98]))
        ws.write(rows, 103, "{}".format(l[99])) #MetaSVM_score
        if l[100] == "D" : #MetaSVM_pred
            ws.write(rows, 104, "{}".format(l[100]), fondo_rojo)
        elif l[100] == "T" :
            ws.write(rows, 104, "{}".format(l[100]), fondo_verde)
        else :
            ws.write(rows, 104, "{}".format(l[100]))
        ws.write(rows, 105, "{}".format(l[101])) #MetaLR_score
        if l[102] == "D" : #MetaLR_pred
            ws.write(rows, 106, "{}".format(l[102]), fondo_rojo)
        elif l[102] == "T" :
            ws.write(rows, 106, "{}".format(l[102]), fondo_verde)
        else :
            ws.write(rows, 106, "{}".format(l[102]))
        ws.write(rows, 107, "{}".format(l[103])) #integrated_fitCons_score
        ws.write(rows, 108, "{}".format(l[104])) #integrated_confidence_value
        ws.write(rows, 109, "{}".format(l[105])) #GERP++_RS
        ws.write(rows, 110, "{}".format(l[106])) #phyloP7way_vertebrate
        ws.write(rows, 111, "{}".format(l[107])) #phyloP20way_mammalian
        ws.write(rows, 112, "{}".format(l[108])) #phastCons7way_vertebrate
        ws.write(rows, 113, "{}".format(l[109])) #phastCons20way_mammalian
        ws.write(rows, 114, "{}".format(l[110])) #SiPhy_29way_logOdds
        ws.write(rows, 115, "{}".format(l[111])) #Filtro variant caller
        ws.write(rows, 116, "{}".format(l[112])) #Columna INFO del variant caller
        ws.write(rows, 117, "{}".format(l[113])) #Columna FORMAT del variant caller
        ws.write(rows, 118, "{}".format(l[114])) #Valores encontrados en la muestra

        if len(l) > 115 :
            ws.write(rows, 119, "{}".format(l[115])) #Valores encontrados en la otra muestra (caso tumor vs normal)
            #Repintar las cabeceras si hay mas de un grupo de valores en la columna format
            ws.write(1, 118, "FORMAT VALUES from VCF(Normal sample)", fondo_azul)
            ws.write(1, 119, "FORMAT VALUES from VCF(Tumor sample)", fondo_azul)

        rows += 1

def crearCabecera(f, wb) :
    '''Crear la cabecera del Excel en la pestana pasada por parametro'''
    #Estilos de la cabecera
    fondo_azul = wb.add_format({'bold' : True,
        'align' : 'center',
        'border' : 1,
        'bg_color' :  '#B3E6FF',
        'font_size' : 13
    })

    #Escribir las cabeceras en cada una de las pestanas
    cols = 0
    for l in superCabecera :
        if l[1] == 1 :
            f.write(0, cols, "{}".format(l[0]), fondo_azul)
        else :
            f.merge_range(0, cols, 0, cols + l[1] -1, "{}".format(l[0]), fondo_azul)
        cols += l[1]

    cols = 0
    for n in cabecera :
        f.write(1, cols, n, fondo_azul)
        cols += 1

def escribirAyuda(nom, i) :
    '''Escribir un cuadro ayuda para interpretar los datos de los predictores de patogenicidad
    @param {String} Nombre del Excel donde se guardara el cuadro
    @param {String} Identificador de las muestras para navegar entre las tablas'''
    wb = openpyxl.load_workbook(filename=nom + ".xlsx")
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
    for it in range(1,5) :
        sh = wb[i+sufijos[it]]
        fila = 2
        print "\tINFO\t{}\tBuscant linia per escriure ajuda".format(getTemps())
        while sh.cell(row=fila, column=8).value != None :
            fila += 1
        print "\tINFO\t{}\tLinia trobada -> {}".format(getTemps(),fila)
        fila += 2 #Margen para escribir el cuadro
        #Escribir la tabla
        sh.cell(row=fila, column=4, value=titol).border = Border(top=Side(style="medium"), left=Side(style="medium"), right=Side(style="medium"))
        sh.cell(row=fila, column=4).font = Font(bold=True)
        sh.merge_cells(start_row=fila, start_column=4, end_row=fila, end_column=9)
        sh.row_dimensions[fila].height = 20
        sh.cell(row=fila+1, column=4, value=sift).border= Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+1, start_column=4, end_row=fila+1, end_column=9)
        sh.row_dimensions[fila+1].height = 40
        sh.cell(row=fila+2, column=4, value=polyphen).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+2, start_column=4, end_row=fila+2, end_column=9)
        sh.row_dimensions[fila+2].height = 50
        sh.cell(row=fila+3, column=4, value=lrt).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+3, start_column=4, end_row=fila+3, end_column=9)
        sh.row_dimensions[fila+3].height = 50
        sh.cell(row=fila+4, column=4, value=mtaster).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+4, start_column=4, end_row=fila+4, end_column=9)
        sh.row_dimensions[fila+4].height = 50
        sh.cell(row=fila+5, column=4, value=massessor).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+5, start_column=4, end_row=fila+5, end_column=9)
        sh.row_dimensions[fila+5].height = 50
        sh.cell(row=fila+6, column=4, value=provean).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+6, start_column=4, end_row=fila+6, end_column=9)
        sh.row_dimensions[fila+6].height = 40
        sh.cell(row=fila+7, column=4, value=fathmm).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+7, start_column=4, end_row=fila+7, end_column=9)
        sh.row_dimensions[fila+7].height = 40
        sh.cell(row=fila+8, column=4, value=msvm).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+8, start_column=4, end_row=fila+8, end_column=9)
        sh.row_dimensions[fila+8].height = 40
        sh.cell(row=fila+9, column=4, value=mlr).border = Border(left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+9, start_column=4, end_row=fila+9, end_column=9)
        sh.row_dimensions[fila+9].height = 40
        sh.cell(row=fila+10, column=4, value=ot).border = Border(bottom=Side(style="medium"), left=Side(style="medium"), right=Side(style="medium"))
        sh.merge_cells(start_row=fila+10, start_column=4, end_row=fila+10, end_column=9)
        sh.row_dimensions[fila+10].height = 40

    wb.save(nom + ".xlsx")

def main() :
    '''Programa principal. Comprueba que el fichero pasado por parametro (sys.argv) sea correcto. Lo lee linea por linea y va llamando a las distintas funciones para trabajar con el archivo'''
    #Parche para intercalar en caso de que la variante de lancet no
    na = -1
    prx = {'CADD_1000g_all' : na, 'CADD_1000g_afr' : na, 'CADD_1000g_amr' : na, 'CADD_1000g_eur' : na, 'CADD_1000g_eas' : na, 'CADD_1000g_sas' : na,
    'dbNSFP_1000g_all' : na, 'dbNSFP_1000g_afr' : na, 'dbNSFP_1000g_amr' : na, 'dbNSFP_1000g_eur' : na, 'dbNSFP_1000g_eas' : na, 'dbNSFP_1000g_sas' : na,
    'CADD_ESP6500_all' : na, 'CADD_ESP6500_ea' : na, 'CADD_ESP6500_aa' : na,
    'dbNSFP_esp6500_all' : na, 'dbNSFP_esp6500_ea' : na, 'dbNSFP_esp6500_aa' : na,
    'ExAC_ExAC_all' : na, 'ExAC_ExAC_afr' : na, 'ExAC_ExAC_amr' : na, 'ExAC_ExAC_eas' : na, 'ExAC_ExAC_fin' : na, 'ExAC_ExAC_nfe' : na, 'ExAC_ExAC_oth' : na, 'ExAC_ExAC_sas' : na,
    'dbNSFP_ExAC_all' : na, 'dbNSFP_ExAC_afr' : na, 'dbNSFP_ExAC_amr' : na, 'dbNSFP_ExAC_eas' : na, 'dbNSFP_ExAC_fin' : na, 'dbNSFP_ExAC_nfe' : na, 'dbNSFP_ExAC_oth' : na, 'dbNSFP_ExAC_sas' : na,
    'dbSNP_MAF' : na}

    if len(sys.argv) < 3 :
        print "\nERROR. No hi ha arxiu per llegir i convertir a Excel.\n\tUs python WES-allInfo2Excel.py variantCall.table_annovar variant_caller(lancet,strelka)\n"
    elif not os.path.isfile(sys.argv[1]) :
        print "\nERROR l'arxiu {} no existeix\n".format(sys.argv[1])
    else :
        aux = sys.argv[1].split(".")
        vc = sys.argv[2]
        nomExcel = ".".join(aux[0:-2])
        ident = nomExcel.split(".")[0]
        inici = time.time()
        perExcel = []
        it = 1
        progres = 0
        percent = 0
        linies = sum(1 for line in open(sys.argv[1]))
        print "INFO\t{}\tBuscant informacio per {} variants".format(getTemps(), linies - 1)
        with open(sys.argv[1], "r") as fi :
            for a in fi :
                if it >= 2 :
                    a = a.strip("\n")
                    c = a.split("\t")
                    try :
                        #if vc == "lancet" :
                        if extraureFiltre(c) == "PASS" :
                            b = getWebInfo.getMyVariant(c[0], c[1], c[3], c[4])
                        else :
                            b = prx
                        #else :
                        #    b = getWebInfo.getMyVariant(c[0],c[1],c[3],c[4])
                    except NameError:
                        b = None

                    final = combinaOrdena(c, b, vc, ident)
                    perExcel.append(final)
                    percent = float(it)/float(linies)*100
                if percent >= progres :
                    print "INFO\t{}\t{} variants buscades. {}% completat".format(getTemps(), it, percent)
                    progres += 10
                it += 1

        print "INFO\t{}\tEscrivint les dades en arxius tsv".format(getTemps())
        #TODO Descomentar
        #escribirTsv(ident, perExcel)
        print "INFO\t{}\tArxius tsv escrits".format(getTemps())
        print "INFO\t{}\tEscrivint Excel".format(getTemps())
        escribirExcel(nomExcel, ident, perExcel)
        print "INFO\t{}\tExcel escrit".format(getTemps())
        t = time.time() - inici
        tt = time.strftime("%H:%M:%S", time.gmtime(t))
        print "INFO\t{}\tPrograma finalitzat. Temps transcorregut: {}".format(getTemps(),tt)

if __name__ == "__main__" :
    main()
