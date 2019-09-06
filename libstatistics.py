#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
MAIN: Functions to extract statistics from the LOH processed data
"""

"""
FUNCTIONS
    PENDING TO BE REMODELED
    - print2Bed
    - print2GGplot
    - extractStatistics
    - getVariations: Obtain TPR, TNR, PPV, NPR, FNR, FPR, FDR, FOR, ACC from 2x2 matrix
    - checkLogR
"""

"""
Libraries
"""
import libcomparison as comp
import libconstants as cts
import math

def countsXtool(regs1, regs2 = None) :
    #Calcular els counts que ha reportat cada eina
    t1 = {"A" : 0, "D" : 0, "L" : 0, "N" : 0}
    if regs2 != None :
        t2 = {"A" : 0, "D" : 0, "L" : 0, "N" : 0}
    for a in cts.chromosomes :
        if a in regs1.keys() :
            for i in regs1[a] :
                t1[i[2]] += 1
        if regs2 != None :
            if a in regs2.keys() :
                for j in regs2[a] :
                    t2[i[2]] += 1
    if regs2 != None :
        return (t1, t2)
    else :
        return t1

def calculateCounts(tab) :
    #Calcular les aberracions que ha reportat cada eina segons les regions en comu
    t1 = {"A" : 0, "D" : 0, "L" : 0, "N" : 0}
    t2 = {"A" : 0, "D" : 0, "L" : 0, "N" : 0}
    for a in cts.aberrations :
        cont = 0
        for i in tab[a].keys() :
            cont += tab[a][i]
        t1[a] = cont
    for a in cts.aberrations :
        cont = 0
        for b in cts.aberrations :
            cont += tab[b][a]
        t2[a] = cont

    return (t1, t2)

#Guardar los datos para que puedan ser usados como un bed en el cual se comparan las regiones de cada una de las herramientas, junto con las regiones en comun.
#Con este bed se compara que las regiones calculadas por getFragments se hayan calculado correctamente
def print2Bed(regs1, prog1, regs2, prog2, regs3) :
    bedfile = "{}_{}_all_regions.bed".format(prog1, prog2)
    with open(bedfile, "w") as fi :
        fi.write("track name=\"{0}\" description=\"{0}\" color=0,128,0 visibility=full db=hg38\n".format(prog1))
        for r in cts.chromosomes :
                if r in regs1.keys() :
                    for k in regs1[r] :
                        fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))

        fi.write("\ntrack name=\"{0}\" description=\"{0}\" color=128,0,0 visibility=full db=hg38\n".format(prog2))
        for r in cts.chromosomes :
            if r in regs2.keys() :
                for k in regs2[r] :
                    fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))

        fi.write("\ntrack name=\"Regions split\" description=\"Splitted regions\" color=0,0,128 visibility=full db=hg38\n")
        for r in cts.chromosomes :
            if r in regs3.keys() :
                for k in regs3[r] :
                    fi.write("chr{}\t{}\t{}\n".format(r, k[0], k[1]))
        print "INFO: Created a bed file with the regions corresponding to {} output, {} output, and the regions in common".format(prog1, prog2)

def printTable(dc, prog1, prog2) :
    #Print a 4x4 table. Calculate the contingency
    txt = "\t{}\n".format("\t".join(cts.aberrations))
    for a in cts.aberrations :
        txt += "{}".format(a)
        for b in cts.aberrations :
            txt += "\t{}".format(dc[a][b])
        txt += "\n"
    table = "{}_{}_4x4.tsv".format(prog1, prog2)
    with open(table, "w") as fi :
        fi.write(txt)
    print "INFO: Created a tab-separated file 4x4 table comparing the output for both programs"

def doContingency(dc, counts1, counts2) :
    allab = 0 #Count the total of regions to get the TN value easily
    sts = {}
    for a in cts.aberrations :
        for b in cts.aberrations :
            allab += dc[a][b]

    for a in cts.aberrations :
        tp = dc[a][a]
        fp = counts1[a] - tp
        fn = counts2[a] - tp
        tn = allab - tp - fp - fn
        sts = getConfusionStatistics(tp, fp, fn, tn)
        print "TP={}".format(tp)
        print "FP={}".format(fp)
        print "FN={}".format(fn)
        print "TN={}".format(tn)
        print sts


#TODO remodelar
#Extraer calculos a partir de una tabla de contingencia
def getConfusionStatistics(tp, fp, fn, tn) :
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
     - F1 score
     - Matthews correlation coefficient (MCC)
    from 2x2 confusion matrix
    """
    stats = {}
    tp = float(tp)
    fp = float(fp)
    fn = float(fn)
    tn = float(tn)
    stats["TPR"] = tp/(tp+fn)
    stats["TNR"] = tn/(tn+fp)
    stats["PPV"] = tp/(tp+fp)
    stats["NPV"] = tn/(tn+fn)
    stats["FNR"] = fn/(fn+tp)
    stats["FPR"] = fp/(fp+tn)
    stats["FDR"] = fp/(fp+tp)
    stats["FOR"] = fn/(fn+tn)
    stats["ACC"] = (tp+tn)/(tp+fp+fn+tn)
    stats["F1"] = (2*tp)/(2*tp + fp + fn)
    stats["MCC"] = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return stats


#TODO remodelar
#Preparar los datos para crear un archivo apto para el script doGGplot.R el cual hara un ggplot del copy number de los outputs de las herramientas
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

#TODO remodelar
#Datos para la tabla 3x3
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

#TODO remodelar
#Preparar datos para una comparacion entre logR
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

if __name__ == "__main__" :
    """
        UNIT TEST
    """
    pr1 = "FACETS"
    pr2 = "ascatngs"
    print "Reading FACETS example"
    fa = comp.convert2region("input_examples/facets_comp_cncf.tsv", pr1)
    print "Reading AscatNGS example"
    s = comp.convert2region("input_examples/TCGA-13-0887-01A-01W.copynumber.caveman.csv", pr2)
    print "Read complete. Getting the fragments"
    regs = comp.getFragments(fa, s)
    print "Got fragments. Checking the copy number"
    #NOTE Comprovar que la eina 1 es FACETS i que la eina 2 es ascatNGS
    dc = comp.doComparison(regs, fa, s)
    print "Copy number done. Preparing some statistics"
    print "1) Counts"
    c1, c2 = calculateCounts(dc)
    print "2) Counts per tool"
    counts1, count2 = countsXtool(fa, s)
    print "3) BED file"
    #print2Bed(fa, pr1, s, pr2, regs)
    print "4) 4x4 comparison table"
    #printTable(dc, pr1, pr2)
    print "5) Contingency table"
    doContingency(dc, c1, c2)
