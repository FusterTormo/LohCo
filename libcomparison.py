#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import math
import libextractfile as ex
import libgetters as ge
import libconstants as cts

"""
MAIN: All functions to do a comparison between two LOH files
"""

"""
FUNCTIONS:
    - convert2region: Reads the output LOH files, and converts the data to REGION format
    - getFragments: Gets a list of regions in common between both files
    - checkCopyNumber: Compare both files. Returns the statistics as needed
"""

"""
REGION FORMAT:
{chr : [start, end, copy-number, total_copy_number, low_copy_number, logR_value]],
ploidy : float,
purity : float,
likelyhood : float}
"""

def convert2region(path, filetype) :
    """Read the output file from LOH program and convert the data to REGION format

    Reads the file located in path passed as parameter and converts the data to REGION format. It is necessary to pass the type of output in order to call the corresponding function to extract the data

    Parameters :
        path (str) : Path to the file where the information is stored
        filetype (str) : Type of program that has output the information passed as parameter. Available values: array, facets, ascat, ascatngs, sequenza, cnacs, titan, purple

    Returns :
        dict : A dict where the key is the chromosome name (chr1, ..., chrX, chrY) and the value is a list of the regions in pairs [start, end]. Additionally values of ploidy, purity and likelyhood.
            Look REGION format
    """
    reg = {}
    filetype = filetype.lower()
    if filetype == "array" :
        reg = ex.extractArray(path)
    elif filetype == "facets" :
        reg = ex.extractFacets(path)
    elif filetype == "ascat" or filetype == "ascatngs" :
        reg = ex.extractAscat(path)
    elif filetype == "sequenza" :
        print "ERROR: Sequenza is not implemented yet"
        sys.exit()
    elif filetype == "cnacs" :
        print "ERROR: CNACS is not implemented yet"
        sys.exit()
    elif filetype == "titan" :
        print "ERROR: Titan is not implemented yet"
        sys.exit()
    elif filetype == "purple" :
        print "ERROR: Purple is not implemented yet"
        sys.exit()
    else :
        raise IOError("ERROR: Type of file not found. Cannot continue")
    return reg

def getFragments(l1, l2) :
    """Split the regions in REGION variables to get the regions with the same size of both

    Split the regions of both lists passed as parameter to get other list with all the regions

    Parameters:
        l1 (dict) : Dictionary, in REGION format, with the output of one program
        l2 (dict) : Dictionary, in REGION format, with the output of the other program

    Returns:
        dict : Using the name of the chromosome as key, and in each chromosome, a list of all the regions with the format [start, end]. The format is similar to REGION format,
                but without the ploidy, purity, and likelyhood keys
    """
    allregs = {}
    aux1 = []
    aux2 = []
    tmp1 = []
    tmp2 = []
    last = 0
    maxChr = {"1" : 248956422, "2" : 242193529, "3" : 198295559,
    "4" : 190214555, "5" : 181538259, "6" : 170805979,
    "7" : 159345973, "8" : 145138636, "9" : 138394717,
    "10" : 133797422, "11" : 135086622, "12" : 133275309,
    "13" : 114364328, "14" : 107043718, "15" : 101991189,
    "16" : 90338345, "17" : 83257441, "18" : 80373285,
    "19" : 58617616, "20" : 64444167, "21" : 46709983,
    "22" : 50818468, "X" : 156040895, "Y" : 999999999999999999999999999999999999999999999}
    iters = 0
    for k in  cts.chromosomes:
        if k in l1.keys() :
            tmp1 = list(l1[k])
        else :
            tmp1 = []

        if k in l2.keys() :
            tmp2 = list(l2[k])
        else :
            tmp2 = []
        if tmp1 == [] and tmp2 == [] :
            tmpregs = [[0, maxChr[k]]]
        else :
            tmpregs = []
            aux1 = []
            aux2 = []
            while len(tmp1) > 0 or len(tmp2) > 0 :
                if len(aux1) == 0 :
                    if len(tmp1) > 0 :
                        aux1 = ge.getRegion(tmp1.pop(0))
                    else :
                        aux1.append(maxChr[k])
                if len(aux2) == 0 :
                    if len(tmp2) > 0 :
                        aux2 = ge.getRegion(tmp2.pop(0))
                    else :
                        aux2.append(maxChr[k])
                reg = []
                # Si els dos array estan buits, la regio comuna es la dels dos valors minims
                if len(tmpregs) == 0 :
                    reg = [0, min(aux1[0], aux2[0])]
                    if min(aux1[0], aux2[0]) == aux1[0] :
                        aux1.pop(0)
                    else :
                        aux2.pop(0)
                else :
                    if min(aux1[0], aux2[0]) == aux1[0] :
                        reg = [last, aux1.pop(0)]
                    else :
                        reg = [last, aux2.pop(0)]
                last = reg[1] + 1
                tmpregs.append(reg)
                iters += 1

            while len(aux1) > 0 or len(aux2) > 0 :
                reg = []
                if len(aux1) == 0 :
                    reg = [last, aux2.pop(0)]
                elif len(aux2) == 0 :
                    reg = [last, aux1.pop(0)]
                else :
                    if min(aux1[0], aux2[0]) == aux1[0] :
                        reg = [last, aux1.pop(0)]
                    else :
                        reg = [last, aux2.pop(0)]
                last = reg[1] + 1
                tmpregs.append(reg)

            if last < maxChr[k] :
                tmpregs.append([last, maxChr[k]])

            del(aux1)
            del(aux2)
        allregs[k] = tmpregs

        del(tmpregs)
    return allregs

def doComparison(regions, t1, t2) :
    """
    Checks the copy number for all the fragments in the region passed as parameter (regions) in each of the 2 tools passed as parameter
    Returns a two-dimension dict with the data like a contingency table
    TODO: end the function's documentation
    """
    tab = {"D" : {"D" : 0, "N" : 0, "L" : 0, "A" : 0}, "N" : {"D" : 0, "N" : 0, "L" : 0, "A" : 0}, "L" : {"D" : 0, "N" : 0, "L" : 0, "A" : 0}, "A" : {"D" : 0, "N" : 0, "L" : 0, "A" : 0}}
    cont = 0
    for chr in cts.chromosomes :
        for r in regions[chr] :
            cont += 1
            c1 = ge.getCopyNumber(r, chr, t1)
            c2 = ge.getCopyNumber(r, chr, t2)
            tab[c1][c2] += 1


    print tab
    print cont


def checkCopyNumber_old(regions, t1, t2, name1 = "tool1", name2 = "tool2") :
    """Checks, for each region in the fragments if the region is called as CNA or CNV.
    Outputs 3x3 contingency table with deletions, amplifications and CNN normal of all the regions.
    Outputs counts for each alteration in each program"""
    #Counters for contingency table
    arDel_faDel = 0
    arDel_faNorm = 0
    arDel_faAmp = 0
    arNorm_faDel = 0
    arNorm_faNorm = 0
    arNorm_faAmp = 0
    arAmp_faDel = 0
    arAmp_faNorm = 0
    arAmp_faAmp = 0
    del_FA = 0
    amp_FA = 0
    norm_FA = 0
    del_ar = 0
    amp_ar = 0
    norm_ar = 0
    cont = 0
    for chr, reg in regions.iteritems() :
        for r in reg :
            c1 = getCopyNumber(r, chr, t1)
            c2 = getCopyNumber(r, chr, t2)
            if c1 == 'D' :
                if c2 == 'D' :
                    arDel_faDel += 1
                    del_FA += 1
                    del_ar += 1
                elif c2 == 'N' :
                    arDel_faNorm += 1
                    norm_FA += 1
                    del_ar += 1
                elif c2 == 'A' :
                    arDel_faAmp += 1
                    amp_FA += 1
                    del_ar += 1
                else :
                    print "ERROR: Valor no reconegut per T2 {}".format(c2)
            elif c1 == 'N' :
                if c2 == 'D' :
                    arNorm_faDel += 1
                    norm_ar += 1
                    del_FA += 1
                elif c2 == 'N' :
                    arNorm_faNorm += 1
                    norm_ar += 1
                    norm_FA += 1
                elif c2 == 'A' :
                    arNorm_faAmp += 1
                    norm_ar += 1
                    amp_FA += 1
                else :
                    print "ERROR: Valor no reconegut per T2 {}".format(c2)
            elif c1 == 'A' :
                if c2 == 'D' :
                    arAmp_faDel += 1
                    amp_ar += 1
                    del_FA += 1
                elif c2 == 'N' :
                    arAmp_faNorm += 1
                    amp_ar += 1
                    norm_FA += 1
                elif c2 == 'A' :
                    arAmp_faAmp += 1
                    amp_ar += 1
                    amp_FA += 1
                else :
                    print "ERROR: Valor no reconegut per T2 {}".format(c2)
            else :
                print "ERROR: Valor no reconegut per T1 {}".format(c1)
            cont += 1

    ct = "{}_{}_counts.tsv".format(name1, name2)
    with open(ct, "w") as fi :
        fi.write("\tDel\tNorm\tAmp\n{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n".format(name1, del_ar, norm_ar, amp_ar, name2, del_FA, norm_FA, amp_FA))
    print "INFO: Counts matrices stored as {}".format(ct)
    extractStatistics(arDel_faDel, arNorm_faDel, arAmp_faDel, arDel_faNorm, arNorm_faNorm, arAmp_faNorm, arDel_faAmp, arNorm_faAmp, arAmp_faAmp, name1, name2)


def print2Bed(regs1, prog1, regs2, prog2, regs3) :
    sr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    bedfile = "{}_{}_all_regions.bed".format(prog1, prog2)
    with open(bedfile, "w") as fi :
        fi.write("track name=\"{0}\" description=\"{0}\" color=0,128,0 visibility=full db=hg38\n".format(prog1))
        for r in sr :
            if r in regs1.keys() :
                for k in regs1[r] :
                    fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))

        fi.write("\ntrack name=\"{0}\" description=\"{0}\" color=128,0,0 visibility=full db=hg38\n".format(prog2))
        for r in sr :
            if r in regs2.keys() :
                for k in regs2[r] :
                    fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))

        fi.write("\ntrack name=\"Regions split\" description=\"Splitted regions\" color=0,0,128 visibility=full db=hg38\n")
        for r in sr :
            if r in regs3.keys() :
                for k in regs3[r] :
                    fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))

    print "INFO: Bed file with all the regions of all programs written as {}".format(bedfile)

def print2GGplot(dc1, dc2, prog1, prog2) :
    sr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    longitud = 5 #Tamano que tiene que tener cada lista para poder sacer el CN necesario para crear la tabla para el ggplot
    txtTCN = "chr\tstart\tend\tcn\ttype\n"
    txtLCN = "chr\tstart\tend\tcn\ttype\n"
    tcn = "{}_{}_tab4ggplot_TCN.tsv".format(prog1, prog2)
    lcn = "{}_{}_tab4ggplot_lcn.tsv".format(prog1, prog2)

    for r in sr :
        if r in dc1.keys() :
            for c in dc1[r] :
                if len(c) < longitud :
                    raise IndexError("ERROR: List from program {} is too short to find CN values".format(prog1))
                else :
                    txtTCN += "{}\t{}\t{}\t{}\t{}\n".format(r, c[0], c[1], c[3], prog1)
                    if c[4] != 'NA' : #NOTE ruling out the regions with lcn == 'NA'
                        txtLCN += "{}\t{}\t{}\t{}\t{}\n".format(r, c[0], c[1], c[4], prog1)
    for r in sr :
        if r in dc2.keys() :
            for c in dc2[r] :
                if len(c) < longitud :
                    raise IndexError("ERROR: List from program {} is too short to find CN values".format(prog2))
                else :
                    txtTCN += "{}\t{}\t{}\t{}\t{}\n".format(r, c[0], c[1], c[3], prog2)
                    if c[4] != 'NA' : #NOTE ruling out the regions with lcn == 'NA'
                        txtLCN += "{}\t{}\t{}\t{}\t{}\n".format(r, c[0], c[1], c[4], prog2)


    with open(tcn, "w") as fi :
        fi.write(txtTCN)
    with open (lcn, "w") as fi :
        fi.write(txtLCN)

    print "INFO: Tables for ggplot written as {} and {}. To create the ggplot".format(tcn, lcn)
    print "WARNING: NA values from lcn have been removed"


if __name__ == "__main__" :
    """
        UNIT TEST
    """
    print "Reading FACETS example"
    fa = convert2region("input_examples/facets_comp_cncf.tsv", "FACETS")
    print "Reading AscatNGS example"
    s = convert2region("input_examples/TCGA-13-0887-01A-01W.copynumber.caveman.csv", "ascatngs")
    print "Read complete. Getting the fragments"
    regs = getFragments(fa, s)
    print "Got fragments. Checking the copy number"
    doComparison(regs, fa, s)
