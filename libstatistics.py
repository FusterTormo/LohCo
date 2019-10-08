#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
MAIN: Functions to extract statistics from the LOH processed data
"""

"""
FUNCTIONS
    - countsXtool : Calculates the number of each alteration reported by each tool
    - calculateCounts : Calculates the number of each alteration in each region in common
    - print2Bed : Stores the regions reported by each tool, and the regions once they are put in common
    - printTable : Prints a 4x4 table with the alterations reported by each tool
    - doContingency : Extracts the 2x2 confusion matrix from the original 4x4 table. Calculates TP, TN, FP, FN, and the confusion matrix statistics
    - getConfusionStatistics : Given the TP, FP, FN, and TN, calculates the sensitivity, specificity, accuracy, F1 score, Matthews correlation coefficient, among others
    - jaccardIndex : Calculates the Jaccard Index from the original 4x4 table
    - doGGplotFiles : Stores the information necessary to do a ggplot with the counts of the copy numbers. Executes the Rscript that creates the plot
    - logRcomp : Stores the information necesary to do a correlation between logR outputs. Executes the Rscript that creates the plot
"""

"""
Libraries
"""
import libcomparison as comp
import libconstants as cts
import libgetters as getlib
import math
import subprocess
import os

def countsXtool(regs1, regs2 = None) :
    """Calculate the number of each alteration reported by each tool

    Counts the number of (A)mplifications, (D)eletions, (L)oss of heterozygosity, and (N)ormal regions reported by each tool. Stores the count information in a dictionary. If no second tool is
    provided, only calculates the counts for one tool. In the other case, two dictionaries are returned

    Parameters :
        regs1 (dict) : Dictionary, in REGION format, with the output of one program
        regs2 (dict) : Dictionary, in REGION format, with the output of the other program. Optional parameter

    Returns :
        dict : If only one parameter is given, one dictionary is returned with the counts of each alteration. In case two parameters are given, two dictionaries will be returned, each with the
        counts of each tool
    """
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
    """Calculate the number of alterations of each tool in the regions in common

    Counts the number of (A)mplifications, (D)eletions, (L)oss of heterozygosity, and (N)ormal reported by each tool in the regions splitted to be in common. Information of each tool is stored in
    a dict for each tool.

    Parameters :
        tab (dict) : Two-dimension dict obtained by executing doComparison function

    Returns :
        t1 (dict) : Sum of rows in tab. That corresponds to the number of each aberration by tool1
        t2 (dict) : Sum of columns in tab parameter. That corresponds to the number of each aberration by tool2
    """
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

def print2Bed(regs1, prog1, regs2, prog2, regs3) :
    """Store the regions by each tool, and the regions splitted to be in common in a bed file

    Gets the regions [List of start-end pairs] reported by tool1 (regs1), and store this coordinates in a bed file. Then gets the regions reported by tool2 and adds this information to
    the bed file. Finally, the regions splitted to be in common between both tools are added to this bed file. Bed file is stored as nameTool1_nameTool2_all_regions.bed. This bed is
    useful to compare if the getFragments execution has been done successfully.

    Parameters :
        regs1 (dict) : Program1's output, in REGION format
        prog1 (str) : Name of program 1. It will be used in the name of the resulting bed file, and in the bed header
        regs2 (dict) : Output of program 2, in REGION format
        prog2 (str) : Name of program 2. It will be used in the name of the resulting bed file, and in the bed header
        regs3 (dict) : Splitted regions in common between both tools. Obtained using getFragments function
    """
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

def printTable(dc, prog1, prog2, print2file = True) :
    """Print 4x4 comparison table.

    Prints a 4x4 table that compares the (D)eletions, (L)oss of heterozygosity, (N)ormal copy number, and (A)mplifications reported by each tool in the regions in common.

    Parameters :
        dc (dict) : Two-dimension dict obtained by executing doComparison function
        prog1 (str) : Name of program1. It will be used as table header
        prog2 (str) : Name of program2. It will be used as table hader
        print2file (bool) : If the information has to be stored in a text file, or be returned as a dict. If print2file is True, no variables will be returned. Default: True

    Returns :
        str : If print2file is False, the file content will be returned
        None : If print2file is True, a text file will be created and nothing is going to be returned
    """
    txt = "{}\t\t{}\n".format(prog1, prog2)
    txt += "\t{}\n".format("\t".join(cts.aberrations))
    for a in cts.aberrations :
        txt += "{}".format(a)
        for b in cts.aberrations :
            txt += "\t{}".format(dc[a][b])
        txt += "\n"
    if print2file :
        table = "{}_{}_4x4.tsv".format(prog1, prog2)
        with open(table, "w") as fi :
            fi.write(txt)
        print "INFO: Created a tab-separated file called {}, with a 4x4 table comparing the output of both programs".format(table)
    else :
        return txt

def doContingency(dc, aber = cts.aberrations) :
    """Extract the information from 4x4 aberrations table to get the confusion matrix stats for the aberrations passed as parameter

    Given the 4x4 aberrations table, calculate the number of true positives (TP), false positives (FP), false negatives (FN), and true negatives (TN) for each aberration. The
    aberrations can be modified by passing a list as parameter. Then the function calls getConfusionStatistics to get TPR, TNR, F1 score, Matthews correlation coefficient, and other statistics.
    All this information is returned.

    Parameters :
        dc (dict) : Two-dimension dict obtained by executing doComparison function
        aber (list) : List of aberrations to do the confusion matrix calculations. Default ['A', 'D', 'L', 'N']

    Returns :
        dict : All the statistics calculated by getConfusionStatistics, plus the TP, FP, FN, TN counts for each of the aberrations passed in 'aber' parameter
    """
    allab = 0 #Count the total of regions to get the TN value easily
    stats = {}
    sts = {}
    counts1, counts2 = calculateCounts(dc)
    for a in aber :
        for b in aber :
            allab += dc[a][b]

    for a in aber :
        tp = dc[a][a]
        fp = counts1[a] - tp
        fn = counts2[a] - tp
        tn = allab - tp - fp - fn
        sts = getConfusionStatistics(tp, fp, fn, tn)
        sts["TP"] = tp
        sts["FP"] = fp
        sts["FN"] = fn
        sts["TN"] = tn
        stats[a] = sts

    return stats

#Extraer calculos a partir de una tabla de contingencia
def getConfusionStatistics(tp, fp, fn, tn) :
    """Calculate all the derivations from a confusion matrix passed as parameter

        From 2x2 confusion matrix, gets
         - True Positive Rate (TPR), sensitivity
         - True Negative Rate (TNR), specificity, selectivity
         - Positive Predictive Value (PPV), precision
         - Negative Predictive Value (NPR)
         - False Negative Rate (FNR), miss rate
         - False Positive Rate (FPR), fall out
         - False Discovery Rate (FDR)
         - False Omission Rate (FOR)
         - Accuracy (ACC)
         - F1 score (F1)
         - Matthews correlation coefficient (MCC)
        and stores all the data in a dict that is returned.

    Parameters :
        tp (int) : Number considered as true positive (TP)
        fp (int) : Number considered as false positive (FP)
        fn (int) : Number considered as false negative (FN)
        tn (int) : Number considered as true negative (TN)

    Returns :
        dict : All calculations in a unique dict. The name of each key corresponds to the abbreviation written before
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

def jaccardIndex(mt, aber = cts.aberrations) :
    """ Calculate the Jaccard Index for each aberration

    Calculates the Jaccard Index in the table given as parameter for the abberations list passed as parameter.

    Parameters :
        mt (dict) :  Two-dimension dict obtained by executing doComparison function. This dict counts the number of bases that corresponds to each region, instead of just the region
        aber (list) : List of aberrations where want to calculate the Jaccard Index. Default ['A', 'D', 'L', 'N']

    Returns :
        dict : Dict with all the calculated Jaccard indexes, one for each aberration passed in aber parameter
    """
    total = 0
    dic = {}
    for a in aber :
        for b in aber :
            total += mt[a][b]

    for a in aber :
        jc = float(mt[a][a]) / float(total)
        dic[a] = jc

    return dic

def doGGplotFiles(reg1, reg2, prog1, prog2, plot = "both") :
    """Create a plot with the total copy number and/or a plot with the minor copy number reported by two tools

    Extracts the information needed from the output of both tools, stores the data in a file with the format needed, and calls the Rscript that creates the plots. Output of both tools has
    to be in REGION format. After the run, a textfile with the data, and a png with the plot will be created

    Parameters :
        reg1 (dict) : Output from tool1 in REGION format
        reg2 (dict) : Output from tool2 in REGION format
        prog1 (str) : Name of the program 1. Will be used as legend in the plot
        prog2 (str) : Name of the program 2. Will be used as legend in the plot
        plot (str) : How many plots create. Options available are :
            - "tcn" to plot only total copy number plot
            - "lcn" to plot only low copy number plot
            - "both" to plot tcn and lcn plots. Default option
    """
    txtTCN = "chr\tstart\tend\tcn\ttype\n"
    txtLCN = "chr\tstart\tend\tcn\ttype\n"
    for a in cts.chromosomes :
        if a in reg1.keys() :
            for b in reg1[a] :
                txtTCN += "{}\t{}\t{}\t{}\t{}\n".format(a, b[0], b[1], b[3], prog1)
                txtLCN += "{}\t{}\t{}\t{}\t{}\n".format(a, b[0], b[1], b[4], prog1)

    for a in cts.chromosomes :
        if a in reg2.keys() :
            for b in reg2[a] :
                txtTCN += "{}\t{}\t{}\t{}\t{}\n".format(a, b[0], b[1], b[3], prog2)
                txtLCN += "{}\t{}\t{}\t{}\t{}\n".format(a, b[0], b[1], b[4], prog2)

    #Write the information in one file, or two depending on plot variable
    if plot == "both" or plot == "tcn" :
        filename = "TCN_ggplot.tsv"
        with open(filename, "w") as fi :
            fi.write(txtTCN)
        print "INFO: Created tab-separated file for ggplot called TCN_ggplot.tsv with tcn copy number information from {} and {}".format(prog1, prog2)
        pythonpath = os.path.dirname(os.path.realpath(__file__))
        rpath = "{}/doGGplot.R {}".format(pythonpath, filename)
        command = "Rscript {} {}_{}_tcn.png".format(rpath, prog1, prog2)
        proc1 = subprocess.Popen(command, shell = True)


    if plot == "both" or plot == "lcn" :
        filename = "lcn_ggplot.tsv"
        with open(filename, "w") as fi :
            fi.write(txtLCN)
        print "INFO: Created tab-separated file for ggplot called lcn_ggplot.tsv with lcn copy number information from {} and {}".format(prog1, prog2)
        pythonpath = os.path.dirname(os.path.realpath(__file__))
        rpath = "{}/doGGplot.R {}".format(pythonpath, filename)
        command = "Rscript {} {}_{}_lcn.png".format(rpath, prog1, prog2)
        proc2 = subprocess.Popen(command, shell = True)

    proc1.communicate()
    proc2.communicate()

def logRcomp(regions, t1, t2, name1 = "tool1", name2 = "tool2") :
    """Create a plot that shows the correlation between logR calculation in both tools.

    Extracts the information needed from the programs output for each of the common regions. Stores this information in a text. Finally calls the Rscript that will create the correlation plot.

    Parameters :
        regions (dict) : Splitted regions to extract the logR from the programs output.
        t1 (dict) : Output from program 1 in REGION format
        t2 (dict) : Output from program 2 in REGION format
        name1 (str) : Name of the tool 1. It will be used in the legend of the plot. Default tool1
        name2 (str) : Name of the tool 2. It will be used in the legend of the plot. Default tool2
    """
    txt = "{}\t{}\n".format(name1, name2)
    fil = "{}_{}_logRcomp.tsv".format(name1, name2)
    for chr, reg in regions.iteritems() :
        for r in reg :
            r1 = getlib.getLogR(r, chr, t1)
            r2 = getlib.getLogR(r, chr, t2)
            if r1 != None and r2 != None :
                txt += "{}\t{}\n".format(r1, r2)
    with open(fil, "w") as fi :
        fi.write(txt)
    print "INFO: Created a file with the logR comparison stored as {}".format(fil)
    pythonpath = os.path.dirname(os.path.realpath(__file__))
    rpath = "{}/compareLogR.R {}".format(pythonpath, fil)
    command = "Rscript {}".format(rpath)
    proc = subprocess.Popen(command, shell = True)
    proc.communicate()

if __name__ == "__main__" :
    """
        UNIT TEST
    """
    pr1 = "FACETS"
    pr2 = "ascatngs"
    print "\n\n\t\tWELCOME TO libstatistics.py UNIT TEST\n\t\t-------------------------------------\n"
    print "Reading FACETS example"
    fa = comp.convert2region("input_examples/facets_comp_cncf.tsv", pr1)
    print "Reading AscatNGS example"
    s = comp.convert2region("input_examples/TCGA-13-0887-01A-01W.copynumber.caveman.csv", pr2)
    print "Read complete. Getting the fragments"
    regs = comp.getFragments(fa, s)
    print "Got fragments. Checking the copy number"
    dc = comp.doComparison(regs, fa, s)
    print "Copy number done. Preparing some statistics"
    print "1) Counts"
    c1, c2 = calculateCounts(dc)
    print "2) Counts per tool"
    counts1, count2 = countsXtool(fa, s)
    print "3) BED file"
    print2Bed(fa, pr1, s, pr2, regs)
    print "4) 4x4 comparison table"
    printTable(dc, pr1, pr2)
    print "5) Contingency table"
    print doContingency(dc)
    print "6) Jaccard index"
    jci = comp.doComparison2(regs, fa, s)
    print jaccardIndex(jci)
