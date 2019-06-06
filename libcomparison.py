#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
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
    lcn = 13
    logR = 4
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
                    reg2 = [int(aux[s]), int(aux[e]), 'N', cnv, aux[lcn].strip(), float(aux[logR])]
                elif cnv > 2 :
                    reg2 = [int(aux[s]), int(aux[e]), 'A', cnv, aux[lcn].strip(), float(aux[logR])]
                else :
                    reg2 = [int(aux[s]), int(aux[e]), 'D', cnv, aux[lcn].strip(), float(aux[logR])]

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
    col_cnv2 = 7
    sc = {}
    sc_cn = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split(",")
            chr = aux[col_c]
            cnv = int(aux[col_cnv])
            reg = [int(aux[col_s]), int(aux[col_e])]
            if cnv == 2 :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'N', cnv, aux[col_cnv2].strip()]
            elif cnv > 2 :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'A', cnv, aux[col_cnv2].strip()]
            else :
                reg2 = [int(aux[col_s]), int(aux[col_e]), 'D', cnv, aux[col_cnv2].strip()]
            if chr in sc.keys() :
                sc[chr].append(reg)
                sc_cn[chr].append(reg2)
            else :
                sc[chr] = [reg]
                sc_cn[chr] = [reg2]

    #Search if copynumber.txt is available to get the logR values
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
        if k in l2.keys() :
            tmp2 = l2[k]
        else :
            tmp2 = []
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

def getLogR(reg, chr, dc) :
    logR = None
    if chr in dc.keys() :
        try :
            for r in dc[chr] :
                if r[0] <= reg[0] and r[1] >= reg[1] :
                    logR = r[5]
                    break
        except IndexError :
            logR = None
    return logR


def checkCopyNumber(regions, t1, t2, name1 = "tool1", name2 = "tool2") :
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


def extractStatistics(del_del, norm_del, amp_del, del_norm, norm_norm, amp_norm, del_amp, norm_amp, amp_amp, name1, name2) :
    mt = "{}_{}_matrix3x3.tsv".format(name1, name2)
    dl = "{}_{}_delVSno.tsv".format(name1, name2)
    nm = "{}_{}_normVSno.tsv".format(name1, name2)
    am = "{}_{}_ampVSno.tsv".format(name1, name2)
    with open(mt, "w") as fi :
        fi.write("\tDel_{0}\tNorm_{0}\tAmp_{0}\n".format(name1))
        fi.write("Del_{}\t{}\t{}\t{}\n".format(name2, del_del, norm_del, amp_del))
        fi.write("Norm_{}\t{}\t{}\t{}\n".format(name2, del_norm, norm_norm, norm_amp))
        fi.write("Del_{}\t{}\t{}\t{}\n".format(name2 ,del_amp, norm_amp, amp_amp))
    print "INFO: 3x3 matrix stored as {}".format(mt)
    with open(dl, "w") as fi :
        fi.write(getVariations(del_del, del_norm+del_amp, norm_del+amp_del, norm_norm+norm_amp+amp_norm+amp_amp))
        fi.write("\n\n{}\t{}\n".format(del_del, del_norm+del_amp))
        fi.write("{}\t{}\n".format(norm_del+amp_del, norm_norm+norm_amp+amp_norm+amp_amp))

    print "INFO: Deletions vs no deletions confusion matrix stored as {}".format(dl)
    with open(nm, "w") as fi :
        fi.write(getVariations(norm_norm, norm_del+norm_amp, del_norm+amp_norm, del_del+del_amp+amp_del+amp_amp))
        fi.write("\n\n{}\t{}\n".format(norm_norm, norm_del+norm_amp))
        fi.write("{}\t{}\n".format(del_norm+amp_norm, del_del+del_amp+amp_del+amp_amp))

    print "INFO: Normal vs no normal confusion matrix stored as {}".format(nm)
    with open(am, "w") as fi :
        fi.write(getVariations(amp_amp, amp_del+amp_norm, del_amp+norm_amp, del_del+del_norm+norm_del+norm_norm))
        fi.write("\n\n{}\t{}\n".format(amp_amp, amp_del+amp_norm))
        fi.write("{}\t{}\n".format(del_amp+norm_amp, del_del+del_norm+norm_del+norm_norm))

    print "INFO: Amplifications vs no amplifications confusion matrix stored as {}".format(am)

def getVariations(tp, fp, fn, tn) :
    """
    Get
     - True Positive Rate (TPR), sensitivity
     - True Negative Rate (TNR), specificity, selectivity
     - Positive Predictive Value (PPV), precision
     - Negative Predictive Value (NPR)
     - False Negative Rate (FNR), miss rate
     - False Positive Rate (FPR), fall out
     - False Discovery Rate (FDR)
     - False Omission Rate (FOR)
     - Accuracy (ACC)
    from 2x2 confusion matrix
    """
    tp = float(tp)
    fp = float(fp)
    fn = float(fn)
    tn = float(tn)
    tpr = tp/(tp+fn)
    tnr = tn/(tn+fp)
    ppv = tp/(tp+fp)
    npv = tn/(tn+fn)
    fnr = fn/(fn+tp)
    fpr = fp/(fp+tn)
    fdr = fp/(fp+tp)
    fomr = fn/(fn+tn)
    acc = (tp+tn)/(tp+fp+fn+tn)
    return "TPR\tTNR\tPPV\tNPV\tFNR\tFPR\tFDR\tFOR\tACC\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(tpr, tnr, ppv, npv, fnr, fpr, fdr, fomr, acc)


def checkLogR(regions, t1, t2, name1 = "tool1", name2 = "tool2") :
    txt = "{}\t{}\n".format(name1, name2)
    fil = "{}_{}_logRcomp.tsv".format(name1, name2)
    for chr, reg in regions.iteritems() :
        for r in reg :
            r1 = getLogR(r, chr, t1)
            r2 = getLogR(r, chr, t2)
            if r1 != None and r2 != None :
                txt += "{}\t{}\n".format(r1, r2)
    with open(fil, "w") as fi :
        fi.write(txt)
    print "INFO: File with the logR comparison stored as {}".format(fil)

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
