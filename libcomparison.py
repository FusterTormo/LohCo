#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import math

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

def getFromFile(path, c, s, e) :
    """Get from the file passed as parameter the chromosome, start, and end position for each line. Return all the data inside a dictionary, using the chromosome as key"""
    chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    ar = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            if chr == "23" :
                chr = "X"
            #Skip headers or columns without interesting information
            if chr in chromosomes :
                reg = [int(aux[s]), int(aux[e])]
                if chr in ar.keys() :
                    ar[chr].append(reg)
                else :
                    ar[chr] = [reg]
    return ar

def getRegionsFromFile(path, c, s, e, v) :
    """Get from the file passed as parameter the chromosome, start, end, and copy_number for each line. Return all the data inside a dictionary, using the chromosome as key"""
    chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    ar = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            if chr == "23" :
                chr = "X"
            #Skip headers or columns without interesting information
            if chr in chromosomes :
                reg = [int(aux[s]), int(aux[e]), float(aux[v])]
                if chr in ar.keys() :
                    ar[chr].append(reg)
                else :
                    ar[chr] = [reg]
    return ar

def extractArray(path) :
    #Get the number of columns where the chromosome, start and end are
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

    #NOTE need to sort the output from each subgroup??
    return (ar, ar_cn)

def extractFacets(path) :
    c = 0
    s = 9
    e = 10
    tcn = 12 #Get tcn column. In case we need lcn column, change to 13
    fa = {}
    fa_cn = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split("\t")
            chr = aux[c]
            if chr == "23" :
                chr = "X"
            #Skip headers or columns without interesting information
            if chr in chromosomes :
                reg = [int(aux[s]), int(aux[e])]
                cnv = int(aux[tcn])
                if cnv == 2 :
                    reg2 = [int(aux[s]), int(aux[e]), 'N', cnv]
                elif cnv > 2 :
                    reg2 = [int(aux[s]), int(aux[e]), 'A', cnv]
                else :
                    reg2 = [int(aux[s]), int(aux[e]), 'D', cnv]

                if chr in fa.keys() :
                    fa[chr].append(reg)
                    fa_cn[chr].append(reg2)
                else :
                    fa[chr] = [reg]
                    fa_cn[chr] = [reg2]

    #NOTE need to sort the output from each subgroup??
    return (fa, fa_cn)

def extractAscat(path) :
    col_c = 1
    col_s = 2
    col_e = 3
    col_cnv = 6 #Get tcn column. In case we need lcn column, change to 7
    sc = {}
    sc_cn = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split(",")
            chr = aux[col_c]
            cnv = int(aux[col_cnv])
            reg = [int(aux[col_s]), int(aux[col_e])]
            if cnv == 2 :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'N', cnv]
            elif cnv > 2 :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'A', cnv]
            else :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'D', cnv]
            if chr in sc.keys() :
                sc[chr].append(reg)
                sc_cn[chr].append(reg2)
            else :
                sc[chr] = [reg]
                sc_cn[chr] = [reg2]

    #NOTE need to sort the output from each subgroup??
    return (sc, sc_cn)

def getFragments(l1, l2) :
    """Split the regions of both lists passed as parameter to get other list with all the regions"""
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
    for k in l1.keys() :
        tmp1 = l1[k]
        tmp2 = l2[k]
        tmpregs = []
        aux1 = []
        aux2 = []
        while len(tmp1) > 0 or len(tmp2) > 0 :
            if len(aux1) == 0 :
                if len(tmp1) > 0 :
                    aux1 = tmp1.pop(0)
                else :
                    aux1.append(maxChr[k])
            if len(aux2) == 0 :
                if len(tmp2) > 0 :
                    aux2 = tmp2.pop(0)
                else :
                    aux2.append(maxChr[k])
            reg = []
            #1) Si els dos array estan buits, la regio comuna es la dels dos valors minims
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

        allregs[k] = tmpregs

        del(tmpregs)
        del(aux1)
        del(aux2)


    return allregs

def getCopyNumber(reg, chr, dc) :
    cnv = "N"
    if chr in dc.keys() :
        for r in dc[chr] :
            if r[0] <= reg[0] and r[1] >= reg[1] :
                cnv = r[2]
                break
    return cnv

def checkCopyNumber(regions, t1, t2, name1 = "tool1", name2 = "tool2") :
    """Checks, for each region in the fragments if the region is called as CNA or CNV"""
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

    print "\nSUMMARY: {} regions analyzed\n\n\t\t3x3 contingency table".format(cont)
    print "\tDel_{0}\tNorm_{0}\tAmp_{0}".format(name1)
    print "Del_{}\t{}\t{}\t{}".format(name2, arDel_faDel, arNorm_faDel, arAmp_faDel)
    print "Norm_{}\t{}\t{}\t{}".format(name2, arDel_faNorm, arNorm_faNorm, arAmp_faNorm)
    print "Del_{}\t{}\t{}\t{}".format(name2 ,arDel_faAmp, arNorm_faAmp, arAmp_faAmp)
    print "\n\t\tCOUNTS\n\tDel\tNorm\tAmp\n{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}".format(name1, del_ar, norm_ar, amp_ar, name2, del_FA, norm_FA, amp_FA)

def print2Bed(regs1, prog1, regs2, prog2, regs3) :
    sr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    bedfile = "all_regions.bed"
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
