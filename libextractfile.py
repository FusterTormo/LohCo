#!/usr/bin/python
# -*- coding: utf-8 -*-

import libconstants as lc

def getFromFile(path, c, s, e) :
    """Read from a file and return the data in a specific format

    Get from the file passed as parameter the chromosome, start, and end position for each line. Return all the data inside a dictionary, using the chromosome as key

    Parameters :
        path (str) : the path of the file to extract the data
        c (int) : column number in the file where the chromosome information is
        s (int) : column number in the file where the start position is
        e (int) : column number in the file where the end position is

    Returns :
        dict : A dict where the key is the chromosome name (chr1, ..., chrX, chrY) and the value is a list of the regions in pairs [start, end]
    """
    ar = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            if chr == "23" :
                chr = "X"
            #Skip headers or columns without interesting information
            if chr in lc.chromosomes :
                reg = [int(aux[s]), int(aux[e])]
                if chr in ar.keys() :
                    ar[chr].append(reg)
                else :
                    ar[chr] = [reg]
    return ar

def extractArray(path) :
    """Read TCGA copy number segment files. Return the data in a specific format

    Read the copy number segment files, extracting the chromsome, start position, end position, and copy number log2 ratio. As these files give less information than other tools, the list is appended with
    'NA' to get the same length as the other tools. So the structure of the list packed in the returning dict will be:
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal
        [3] 'NA' as here is expected the tcn value
        [4] 'NA' as here is expected the lcn value
        [5] 'NA' as here is expected the logR value

    Parameters:
        path (str) : Path of the file to extract the data

    Returns :
        dict : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, 'NA', 'NA', 'NA']
    """
    #Get the number of columns where the chromosome, start and end are
    #TODO Agregar las columnas que le faltan a la lista
    #TODO Eliminar la variable ar_cn
    c = 1
    s = 2
    e = 3
    cnv = 5
    ar = {}
    ar_cn = {}
    #Open the file and
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            #Skip headers or columns without interesting information
            if chr in chromosomes :
                reg = [int(aux[s]), int(aux[e])]
                va = float(aux[cnv])
                if va == 0 :
                    reg2 = [int(aux[s]), int(aux[e]), 'N']
                elif va > 0 :
                    reg2 = [int(aux[s]), int(aux[e]), 'A']
                else :
                    reg2 = [int(aux[s]), int(aux[e]), 'D']
                if chr in ar.keys() :
                    ar[chr].append(reg)
                    ar_cn[chr].append(reg2)
                else :
                    ar[chr] = [reg]
                    ar_cn[chr] = [reg2]

    return (ar, ar_cn)

def extractFacets(path) :
    """Read FACETS $cncf file table and return the interesting information in a specific variable format

    Read FACETS *_cncf.tsv file, which stores the raw information about copy number that has been calculated. From this file the information is converted to a list of regions for each chromosome. The
    structure of the list is as this
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal, 'L' if there is
            a copy neutral LOH
        [3] total copy number (tcn) (int)
        [4] low copy number (lcn) (int)
        [5] logR (float) logR calculation by FACETS
    If tcn is bigger that 2, the region is considered (A)mplified
    If tcn==2, and lcn==1, the region is considered (N)ormal
    If tcn==2, and lcn==0, the region is considered Copy Number Neutral (L)oss of Heterozygosity
    If tcn<2, and lcn<=1, the region is considered (D)eleted

    Parameters :
        path (str) : Path of the file to extract the data

    Returns :
        dict : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]
    """
    c = 0 #Chromosome column
    s = 9 #Start column
    e = 10 #End column
    t = 12 #tcn column
    l = 13 #lcn column
    lR = 4 #logR column
    fa = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            tcn = int(aux[t])
            lcn = int(aux[l])
            logR = float(aux[lR])
            if chr == "23" :
                chr = "X"
            #Skip headers or columns without interesting information
            if chr in lc.chromosomes :
                reg = [int(aux[s]), int(aux[e]), getCN(tcn, lcn), tcn, lcn, logR]

                if chr in fa.keys() :
                    fa[chr].append(reg)
                else :
                    fa[chr] = [reg]
            else :
                print "WARNING: Chromosome {} not found in the chromosomes constant".format(chr)
    return fa


def extractAscat(path) :
    """Read ascatNGS data and return the information in a specific format

    Read ascatNGS *copynumber.caveman.txt file, which stores the raw information about copy number that has been calculated. From this file the information is converted a list of regions for each chromosome. The
    structure of the list is as this
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal, 'L' if there is
            a copy neutral LOH
        [3] total copy number (tcn) (int)
        [4] low copy number (lcn) (int)
        [5] logR (float) logR calculation by FACETS
    If tcn is bigger that 2, the region is considered (A)mplified
    If tcn==2, and lcn==1, the region is considered (N)ormal
    If tcn==2, and lcn==0, the region is considered Copy Number Neutral (L)oss of Heterozygosity
    If tcn<2, and lcn<=1, the region is considered (D)eleted
    LogR is extracted from *copynumber.txt file. Function checks if the file exists. In case the file does not exist, then the logR is completed with 'NA' string

    Parameters :
        path (str) : Path of the file to extract the data

    Returns :
        dict : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]
    """
    #TODO Refer tota la funcio per adaptar al nou format de les variables
    #TODO Eliminar la variable sc_cn
    col_c = 1
    col_s = 2
    col_e = 3
    col_tcn = 6 #Get tcn column. In case we need lcn column, change to 7
    col_lcn = 7
    sc = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split(",")
            chr = aux[col_c]
            tcn = int(aux[col_tcn])
            lcn = int(aux[col_lcn])
            reg = [int(aux[col_s]), int(aux[col_e]), getCN(tcn, lcn), tcn, lcn]

            if chr in sc.keys() :
                sc[chr].append(reg)
            else :
                sc[chr] = [reg]

    #Search if copynumber.txt is available to get the logR values
    #TODO continuar per aci per extraure la informacio de logR des de copynumer.txt
    path2 = path.split(".copynumber")[0] + ".copynumber.txt"
    cr = 1
    pos = 2
    logR = 4
    mtLogR = {}
    cab = False
    if os.path.isfile(path2) :
        print "INFO: Extracting info from {}".format(path2)
        with open(path2, "r") as fi :
            for l in fi :
                if not cab :
                    cab = True
                else :
                    aux = l.split("\t")
                    crom = aux[cr]
                    if crom == "23" :
                        crom = "X"
                    if aux[cr] in mtLogR.keys() :
                        mtLogR[crom][int(aux[pos])] = float(aux[logR])
                    else :
                        mtLogR[crom] = {int(aux[pos]) : float(aux[logR])}
        for cr in sc_cn :
            for ps in sc_cn[cr] :
                if cr in mtLogR.keys() :
                    if ps[0] in mtLogR[cr].keys() :
                        ps.append(mtLogR[cr][ps[0]])
                    elif ps[1] in mtLogR[cr].keys() :
                        ps.append(mtLogR[cr][ps[1]])
                    else :
                        print "WARNING: Region {}:{}-{} not found in {}".format(cr, ps[0], ps[1], path2)

    else :
        print "WARNING: {} not found. LogR calculations could not be added to ASCAT".format(path2)

    #NOTE need to sort the output from each subgroup??
    return (sc, sc_cn)

def getCN(tcn, lcn) :
    """Get the copy number aberration, given the total copy number and the low copy number

    Depending on the value of tcn and lcn, returns if there is an (A)mplification, (D)eletion, Copy number (N)ormal, or Copy Number Neutral (L)oss of Heterozygosity

    Parameters :
    tcn (int) : Total copy number of the region
    lcn (int) : Low copy number of the region

    Returns :
    char : 'A' if there is an (Amplification), 'D' if there is a deletion, 'N' if the region is normal, 'L' if there is CNN-LOH
    """
    ret = ""
    if tcn > 2 :
        ret = "A"
    elif tcn == 2 and lcn == 1 :
        ret = "N"
    elif tcn == 2 and lcn == 0 :
        ret = "L"
    elif tcn == 2 and lcn ==  2 :
        ret = "L"
    elif tcn < 2 :
        ret = "D"
        print "INFO: For tcn={}, and lcn={} returned {}".format(tcn, lcn, ret)
        return ret
