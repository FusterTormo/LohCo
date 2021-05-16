#!/usr/bin/python
# -*- coding: utf-8 -*-

import libconstants as lc
import os
import csv

"""
MAIN: All functions that extract information from the output files from each program and converts the data to the REGION format
"""

"""
FUNCTIONS
    - getFromFile: Read from a text file and return the data in region format
    - extractArray: Read TCGA copy number segments and convert the data to region format
    - extractFacets:  Read FACETS file information and convert the data to region format
    - extractAscat: Read ascatNGS file information and convert the data to region format
    - getCN: Get the copy number aberration, given the total copy number and the low copy number
    - getAscatLogR: Check in the ascatNGS copynumber file if there is a SNP in the region passed as parameter
"""

"""
REGION FORMAT:
{chr : [start, end, copy-number, total_copy_number, low_copy_number, logR_value]],
ploidy : float,
purity : float,
likelyhood : float}
"""

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
    """Read TCGA copy number segment files. Return the data in a REGION format

    Read the copy number segment files, extracting the chromosome, start position, end position, and copy number log2 ratio. As these files give less information than other tools, the list is appended with
    'NA' to get the same length as the other tools. So the structure of the list packed in the returning dict will be:
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal
        [3] 'NA' as here is expected the tcn value
        [4] 'NA' as here is expected the lcn value
        [5] 'NA' as here is expected the logR value
    The additional information, like purity, is not available. So this information is given as "NA", too.

    Parameters:
        path (str) : Path of the file to extract the data

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, 'NA', 'NA.', 'NA']
    """
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
            if chr in lc.chromosomes :
                va = float(aux[cnv])
                if va == 0 :
                    reg = [int(aux[s]), int(aux[e]), 'N', 'NA', 'NA.', 'NA'] #Array regions does not have information about total_copy_number, low_copy_number, and logR
                elif va > 0 :
                    reg = [int(aux[s]), int(aux[e]), 'A', 'NA', 'NA.', 'NA'] #Array regions does not have information about total_copy_number, low_copy_number, and logR
                else :
                    reg = [int(aux[s]), int(aux[e]), 'D', 'NA', 'NA.', 'NA'] #Array regions does not have information about total_copy_number, low_copy_number, and logR
                if chr in ar.keys() :
                    ar[chr].append(reg)
                else :
                    ar[chr] = [reg]

    #Extra information is not available from TCGA array files
    ar["purity"] = 'NA'
    ar["ploidy"] = 'NA'
    ar["likelyhood"] = 'NA'
    return ar

def extractAscatArray(path) :
    """Read TCGA allele-specific copy number segment files. Return the data in a REGION format

    Read the allele-specific copy number segment files, extracting the chromosome, start position, end position, total copy number, and minor
    copy number data. As these files give less information than other tools, the list is appended with 'NA' to get the same length as the other
    tools. So the structure of the list packed in the returning dict will be:
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal
        [3] total copy number, extracted from the Copy_Number column (int)
        [4] low copy number, extracted from the Minor_Copy_Number column (int)
        [5] 'NA' as here is expected the logR value
    The additional information, like purity, is not available. So this information is given as "NA", too.

    Parameters:
        path (str) : Path of the file to extract the data

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, 'NA', 'NA.', 'NA']
    """
    #Get the number of columns where the chromosome, start and end are
    c = 1 # chromosome column number
    s = 2 # start column number
    e = 3 # end column number
    t = 4 # total copy number column number
    l = 6 # low copy number column number
    ar = {}
    tcn = -1
    lcn = -1
    chr = ""
    # Open the file and convert the data to REGION format
    with open(path, "r") as fi :
        for it in fi :
            aux = it.strip().split("\t")
            chr = aux[c].replace("chr", "")
            # Skip the header
            if chr in lc.chromosomes :
                tcn = int(aux[t])
                lcn = int(aux[l])
                reg = [int(aux[s]), int(aux[e]), getCN(tcn, lcn), tcn, lcn, 'NA']
                if chr in ar.keys() :
                    ar[chr].append(reg)
                else :
                    ar[chr] = [reg]
    ar["likelyhood"] = 'NA'
    ar["purity"] = 'NA'
    ar["ploidy"] = 'NA'
    return ar

def extractFacets(path, verbosity = "warning") :
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
    Additionally, it checks if *_basic.tsv file exists, to get the additional information to the REGION variable

    Parameters :
        path (str) : Path of the file to extract the data
        verbosity (str) : If the function must be verbose: "warning" means that warning messages must be shown. Otherwise, not

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]. Additionally, "ploidy", "purity",
            and "likelyhood" information in separated keys.
    """
    c = 0 #Chromosome column
    s = 9 #Start column
    e = 10 #End column
    t = 12 #tcn column
    l = 13 #lcn column
    lR = 4 #logR column
    fa = {}
    with open(path, "r") as fi :
        for lin in fi :
            if not lin.strip("\"").startswith("chrom") :
                aux = lin.strip("\n").split("\t")
                chr = aux[c]
                tcn = int(aux[t])
                if aux[l] == "NA" :
                    lcn = -1
                else :
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
                elif verbosity == "warning" :
                    print("WARNING: Chromosome {} not found in the chromosomes constant".format(chr))

    fa["likelyhood"] = 'NA'
    fa["purity"] = 'NA'
    fa["ploidy"] = 'NA'
    basicPath = path.replace("_cncf.tsv", "_basic.tsv")
    if os.path.isfile(basicPath) :
        with open(basicPath, "r") as fi :
            header = fi.readline()
            body = fi.readline()
        data = body.split("\t")
        fa["likelyhood"] = float(data[0])
        fa["purity"] = float(data[1])
        fa["ploidy"] = float(data[2])

    if fa["ploidy"] == 'NA' and verbosity == "warning":
        print("WARNING: {} not found. Information about ploidy, purity, and likelyhood not given".format(basicPath))

    return fa

def extractAscat(path, verbosity = "warning") :
    """Read ascatNGS data and return the information in a specific format

    Read ascatNGS *copynumber.caveman.txt file, which stores the raw information about copy number that has been calculated. From this file the information is converted a list of regions for each chromosome. The
    structure of the list is as this
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal, 'L' if there is
            a copy neutral LOH
        [3] total copy number (tcn) (int)
        [4] low copy number (lcn) (int)
        [5] logR (float) logR calculation
    If tcn is bigger that 2, the region is considered (A)mplified
    If tcn==2, and lcn==1, the region is considered (N)ormal
    If tcn==2, and lcn==0, the region is considered Copy Number Neutral (L)oss of Heterozygosity
    If tcn<2, and lcn<=1, the region is considered (D)eleted
    LogR median is extracted from *copynumber.txt file. Function checks if the file exists. In case the file does not exist, then the logR is completed with 'NA' string
    Additionally, it checks if *.samplestatistics.txt file exists, to get the additional information to the REGION variable

    Parameters :
        path (str) : Path of the file to extract the data
        verbosity (str) : If the function must be verbose: "warning" means that warning messages must be shown. Otherwise, not

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]. Additionally, "ploidy", "purity",
            and "likelyhood" information in separated keys.
    """
    col_c = 1
    col_s = 2
    col_e = 3
    col_tcn = 6 #Get tcn column. In case we need lcn column, change to 7
    col_lcn = 7
    sc = {}
    mtLogR = {}
    with open(path, "r") as fi :
        for l in fi :
            aux = l.split(",")
            chr = aux[col_c]
            tcn = int(aux[col_tcn])
            lcn = int(aux[col_lcn])
            reg = [int(float(aux[col_s])), int(float(aux[col_e])), getCN(tcn, lcn), tcn, lcn]

            if chr in sc.keys() :
                sc[chr].append(reg)
            else :
                sc[chr] = [reg]

    #Search if copynumber.txt is available to get the logR values
    path2 = path.replace(".copynumber.caveman.csv", ".copynumber.txt")
    cr = 1
    pos = 2
    logR = 4

    cab = False
    if os.path.isfile(path2) :
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

            for cr in sc :
                for ps in sc[cr] :
                    if cr in mtLogR.keys() :
                        if ps[0] in mtLogR[cr].keys() :
                            ps.append(mtLogR[cr][ps[0]])
                        elif ps[1] in mtLogR[cr].keys() :
                            ps.append(mtLogR[cr][ps[1]])
                        else :
                            pass
                            #NOTE This function slows down a lot the execution. If possible avoid it
                            #getAscatLogR(path2, ps, cr)
    elif verbosity == "warning" :
        print("WARNING: {} not found. LogR calculations could not be added to ASCAT".format(path2))

    path3 = path.replace(".copynumber.caveman.csv", ".samplestatistics.txt")
    sc["likelyhood"] = 'NA'
    sc["purity"] = 'NA'
    sc["ploidy"] = 'NA'
    if os.path.isfile(path3) :
        with open(path3, "r") as fi :
            for l in fi :
                aux = l.split(" ")
                if aux[0].startswith("NormalContamination") :
                    sc["purity"] = 1 - float(aux[1])
                if aux[0].startswith("Ploidy") :
                    sc["ploidy"] = float(aux[1])
                if aux[0].startswith("goodnessOfFit") :
                    sc["likelyhood"] = float(aux[1])

    elif verbosity == "warning" :
        print("WARNING: {} not found. Ploidy, purity, and goodness of fit data could not be added to ASCAT".format(path3))
    return sc

def extractSequenza(path, verbosity = "warning") :
    """Read Sequenza data and return the information in a specific format

    Read Sequenza *_segments.txt file, which stores the raw information about copy number that has been calculated. From this file the information is converted to a list of regions for each chromosome. The
    structure of the list is as this
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal, 'L' if there is
            a copy neutral LOH
        [3] total copy number (tcn) (int). In Sequenza, this value is represented by the CNt column
        [4] low copy number (lcn) (int). In Sequenza, this value is represented by the B column
        [5] NA as Sequenza does not report logR
    If tcn is bigger that 2, the region is considered (A)mplified
    If tcn==2, and lcn==1, the region is considered (N)ormal
    If tcn==2, and lcn==0, the region is considered Copy Number Neutral (L)oss of Heterozygosity
    If tcn<2, and lcn<=1, the region is considered (D)eleted
    Additionally, it checks if *.confints_CP.txt file exists, to get the additional information to the REGION variable

    Parameters :
        path (str) : Path of the file to extract the data
        verbosity (str) : If the function must be verbose: "warning" means that warning messages must be shown. Otherwise, not

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]. Additionally, "ploidy", "purity",
            and "likelyhood" information in separated keys.
    """
    col_c = 0
    col_s = 1
    col_e = 2
    col_tcn = 9
    col_lcn = 11
    seq = {}
    with open(path, "r") as fi :
        for l in fi :
            if not l.startswith("chromosome") :
                aux = l.strip("\n").split("\t")
                chr = aux[col_c].replace("chr", "")
                if aux[col_tcn] == "NA" :
                    tcn = -1
                else :
                    tcn = int(aux[col_tcn])
                if aux[col_lcn] == "NA" :
                    lcn = -1
                else :
                    lcn = int(aux[col_lcn])
                if chr in lc.chromosomes :
                    reg = [int(float(aux[col_s])), int(float(aux[col_e])), getCN(tcn, lcn), tcn, lcn, "NA"]
                    if chr in seq.keys() :
                        seq[chr].append(reg)
                    else :
                        seq[chr] = [reg]
                elif verbosity == "warning" :
                    print("WARNING: Chromosome {} not found in the chromosomes constant".format(chr))

    path2 = path.replace("_segments.txt", "_alternative_solutions.txt")
    seq["ploidy"] = "NA"
    seq["purity"] = "NA"
    seq["likelyhood"] = "NA"
    if os.path.isfile(path2) :
        with open(path2, "r") as fi :
            for l in fi :
                if not l.startswith("cellularity") :
                    aux = l.strip("\n").split("\t")
                    seq["ploidy"] = float(aux[1])
                    seq["purity"] = float(aux[0]) #NOTE: Sequenza uses the term "cellularity", instead of ploidy.
    elif verbosity == "warning" :
        print("WARNING: File {} not found. Cannot add the ploidy and purity values to sequenza variable".format(path2))
    path3 = path.replace("_segments.txt", "_confints_CP.txt")
    if os.path.isfile(path3) :
        with open(path3, "r") as fi :
            for l in fi :
                if not l.startswith("cellularity") :
                    aux = l.strip("\n").split("\t")
                    if float(aux[0]) == seq["purity"] and float(aux[1]) == seq["ploidy"] :
                        seq["likelyhood"] = float(aux[2])

        if seq["likelyhood"] == "NA" and verbosity == "warning":
            print("WARNING: Confidence interval not found in {path} for ploidy {ploidy} and purity {purity}".format(path = path3, ploidy = seq["ploidy"], purity = seq["purity"]))
    elif verbosity == "warning" :
        print("WARNING: File {} not found. Cannot add likelyhood value to sequenza variable".format(path3))

    return seq

def extractPurple(path, verbosity = "warning") :
    """Read PURPLE data and return the information in a specific format

    Read PURPLE *.purple.cnv.somatic.tsv file, which stores the raw information about copy number that has been calculated. From this file the information is converted to a list of regions for each chromosome. The
    structure of the list is as this
        [0] start position of the region (int)
        [1] end position of the region (int)
        [2] copy-number (enum) 'A' if region an amplification is considered in the region, 'D' if a deletion is considered in the region, 'N' if the region is considered copy number normal, 'L' if there is
            a copy neutral LOH
        [3] total copy number (tcn) (int). In PURPLE, this value is represented by the copyNumber column
        [4] low copy number (lcn) (int). In PURPLE, this value is represented by the minorAlleleCopyNumber column
        [5] NA as PURPLE does not report logR
    If tcn is bigger that 2, the region is considered (A)mplified
    If tcn==2, and lcn==1, the region is considered (N)ormal
    If tcn==2, and lcn==0, the region is considered Copy Number Neutral (L)oss of Heterozygosity
    If tcn<2, and lcn<=1, the region is considered (D)eleted
    Additionally, it reads from *.purity.tsv the purity, ploidy and likelyhood data to add it to the REGION variable

    Parameters :
        path (str) : Path of the file to extract the data
        verbosity (str) : If the function must be verbose: "warning" means that warning messages must be shown. Otherwise, not

    Returns :
        REGION : A dict where the key is the chromosome name and the value is a list of the regions in format [start, end, copy-number, tcn, lcn, logR]. Additionally, "ploidy", "purity",
            and "likelyhood" information in separated keys.
    """
    col_c = 0
    col_s = 1
    col_e = 2
    col_tcn = 3
    col_lcn = 14
    pur = {}
    with open(path, "r") as fi :
        for l in fi :
            if not l.startswith("chromosome") :
                aux = l.strip("\n").split("\t")
                chr = aux[col_c].replace("chr", "") # Remove chr prefix
                tcn = round(float(aux[col_tcn])) # Convert copy number to int as it is output as float
                try :
                    lcn = round(float(aux[col_lcn]))
                except :
                    print(aux)
                    sys.exit()
                if chr in lc.chromosomes :
                    reg = [int(float(aux[col_s])), int(float(aux[col_e])), getCN(tcn, lcn), tcn, lcn, "NA"]
                    if chr in pur.keys() :
                        pur[chr].append(reg)
                    else :
                        pur[chr] = [reg]
                elif verbosity == "warning" :
                    print("WARNING: Chromosome {} not found in the chromosomes constant".format(chr))
    path2 = path.replace(".cnv.somatic.tsv", ".purity.tsv")
    pur["purity"] = "NA"
    pur["ploidy"] = "NA"
    pur["likelyhood"] = "NA"
    if os.path.isfile(path2) :
        with open(path2, "r") as fi :
            for l in fi :
                aux = l.strip().split("\t")
                if l.startswith("purity") :
                    puri = aux.index("purity")
                    plo = aux.index("ploidy")
                    lik = aux.index("score")
                else :
                    pur["purity"] = float(aux[puri])
                    pur["ploidy"] = float(aux[plo])
                    pur["likelyhood"] = float(aux[lik])
    elif verbosity == "warning" :
        print("WARNING: File {} not found. Cannot add information regarding purity, ploidy and likelyhood to PURPLE variable".format(path2))

    return pur


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
    elif lcn == -1 : #In case low copy number value is calculated as "NA" by programs like FACETS, we only use the total copy number as a clue
        if tcn == 2 :
            ret = "N"
    else :
        print("WARNING: not found value for tcn -> {} and lcn -> {}".format(tcn, lcn))

    return ret

def getAscatLogR(path, reg, chr) :
    """
        Check in the ascatNGS copynumber file if there is a SNP in the region passed as parameter. This function is used only in case the extractAscat function cannot find the region in the file.
        Please, don't use it as it takes a long of execution time!
    """
    logR = 'NA'
    with open(path) as fi :
        reader = csv.DictReader(fi, delimiter="\t")
        for r in reader :
            auxCr = r['Chromosome']
            if r['Chromosome'] == '23' :
                auxCr = 'X'
            if auxCr == chr and int(r['Position']) >= reg[0] and int(r['Position']) <= reg[1] :
                logR = r['segmented LogR']
                break
    if logR == 'NA' :
        print("Not found region")
    return logR

if __name__ == "__main__" :
    print("INFO:  This is a library file. Don't execute it alone")
