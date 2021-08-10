#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import subprocess
import sys
import time

import main1 as lib
import main8 as asc
import libconstants as cte

"""MAIN PROGRAM
Search VUS variants in LOH
"""

"""
USAGE: Modify the constants below:
    * gene: gene coordinates. Some examples available below
    * genename: gene name to search in ANNOVAR output data
    * varCallSuffix: variant caller of interest. Available options are commented

Additionally, it is possible to modify the TCGA cancer repository by changing the constant cancer (cancer repository name) and cFolder (cancer folder where all the data is)
"""

# Constants
# Gene of interest to search LOH
brca1 = ["17", 43044295, 43125483]
brca2 = ["13", 32315508, 32400268]
atm = ["11", 108222832, 108369099]
palb2 = ["16", 23603165, 23641310]
rad51c = ["17", 58692602, 58735611]
rad51 = ["15", 40695174, 40732340]
# Cancer repository of interest, and full path to files
cFolder = "fsupek_cancer2" # Shared folder where the data is stored
cancer = "OV" # Cancer repository to search
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/{cancer_path}/TCGA_bam/{c}".format(c = cancer, cancer_path = cFolder) # Full path to TCGA repository data
# Constants for ClinVar
keyword = "cancer" # If the variant has this keyword in ClinVar annotations. It will be considered pathogenic
vcf = "clinvar.vcf" # Path to ClinVar database

# Input data. Main file read this variables to do the analysis
gene = palb2
genename = "PALB2"
varCallSuffix = "platypusGerm/platypus.hg38_multianno.txt"
#varCallSuffix = "strelkaGerm/results/variants/strelka.hg38_multianno.txt"

def getTime() :
    """Returns current time (hours:minutes:seconds) in a fancy format"""
    return time.strftime("%H:%M:%S", time.gmtime())

def readClinVar() :
    """Read ClinVar vcf database. Convert the variants and its annotations in a python dict, where the key is the pair chr-position"""
    cnt = {}
    with open(vcf, "r") as fi :
        for l in fi :
            if not l.startswith("#") :
                aux = l.strip().split("\t")
                idx = "{}-{}".format(aux[0], aux[1])
                if idx not in cnt.keys() : # PATCH!! Get the first element if there are position duplicates
                    cnt[idx] = l
    return cnt

def extractClinVar(clinvar, variant) :
    """Extract the important information from ClinVar

    Columns extracted (description given in clinvar vcf file)
        * CLNDISDB: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
        * CLNDN: ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
        * CLNSIG: Clinical significance for this single variant
        * CLNREVSTAT: ClinVar review status for the Variation ID
    """

    data = {"db" : "NA", "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
    # Check if the variant is exactly the same as the reported by ClinVar
    ref = variant[3]
    alt = variant[4]
    aux = variant[0:5]
    # Extract the interesting data from the clinvar line
    tmp = clinvar.split("\t")
    cln_ref = tmp[3]
    cln_alt = tmp[4]
    info = tmp[7]
    aux2 = tmp[0:5]
    # Indel scenario. VCF syntax is different to ANNOVAR syntax. Converting the VCF syntax to ANNOVAR syntax
    if len(cln_ref) > 1 or len(cln_alt) > 1 :
        if len(cln_ref) == 1 :
            cln_ref = "-"
        else :
            cln_ref = cln_ref[1:]
        if len(cln_alt) == 1 :
            cln_alt = "-"
        else :
            cln_alt = cln_alt[1:]

    if ref == cln_ref and alt == cln_alt :
        for i in info.split(";") :
            tmp = i.split("=")
            # Database name identifiers in pairs "OMIM:MMMMMMM"
            if tmp[0] == "CLNDISDB" :
                data["db"] = ",".join(aux[1:])
            # Diagnostic name
            if tmp[0] == "CLNDN" :
                data["disease"] = tmp[1]
            # Clinical significance
            if tmp[0] == "CLNSIG" :
                data["significance"] = tmp[1]
            # Revision status
            if tmp[0] == "CLNREVSTAT" :
                data["revStatus"] = tmp[1]

    return data

def annotateClinVar(data, cln) :
    """Add ClinVar information to the variants list passed as parameter (data)"""
    d = data.split("\t")
    if d[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] :
        anno = []
        chrom = d[0].replace("chr", "")
        # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
        if d[1] == "-" :
            pos = int(d[1]) - 1
        else :
            pos = int(d[1])

        search = "{}-{}".format(chrom, pos)
        supData = {"db" : {}, "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
        if search in cln.keys() :
            supData = extractClinVar(cln[search], d)

        anno = {"chrom" : d[0], "start" : d[1], "end" : d[2], "ref" : d[3], "alt" : d[4], "type" : d[5], "exonicType" : d[8], "disease" : supData["disease"], "significance" : supData["significance"]}
    else :
        anno = ""

    return anno

def saveHistogram(data, filename) :
    """Save dict with number of variants in each gene position in a tsv file"""
    minim = min(data.keys())
    maxim = max(data.keys())
    total = sum(data.values())
    with open(filename, "w") as fi :
        fi.write("position\ttimes\tfreq\n")
        for i in range(minim, maxim+1) :
            if i in data.keys() :
                fi.write("{}\t{}\t{}\n".format(i, data[i], round(data[i]/total, 6)))
            else :
                fi.write("{}\t0\t0\n".format(i))

    print("INFO: Variant histogram data saved as {}".format(filename))

def saveVariants(data, filename) :
    """Save variants data in a tsv file"""
    with open(filename, "w") as fi :
        fi.write("chrom\tstart\tend\tref\talt\ttype\texonic.type\tlohs\tcln.signf\tcln.disease\tsubmitter\ttumor.id\tcontrol.id\n")
        for d in data :
            fi.write("{cr}\t{st}\t{nd}\t{ref}\t{alt}\t{tp}\t{ex}\t{loh}\t{sgnf}\t{dss}\t{sub}\t{tm}\t{cn}\n".format(
                cr = d["chrom"], st = d["start"], nd = d["end"], ref = d["ref"], alt = d["alt"], tp = d["type"], ex = d["exonicType"], loh = d["lohCount"], sgnf = d["significance"], dss = d["disease"],
                sub = d["submitter"], tm = d["tumor"], cn = d["control"]))

    print("INFO: Variants saved as {}".format(filename))

def getData() :
    """Get LOH and variants information for each submitter in the gene of interest

    Classify the submitters as
        * Positive: When 2 or more LOH tools report LOH (CNN-LOH or CNL-LOH) in the gene of interest
        * Pathogenic: Positive submitters that have a pathogenic variant (see cte.var_positive) in combination with LOH
        * Negative: When less than 2 LOH tools report LOH in the gene of interest

    Variants are classified in each list (posData or negData) according to LOH found
    """
    # Variables
    pairs = 0 # Number of submitters with pair tumor-control
    done = [] # List of submitters with ASCAT2, FACETS and Sequenza done
    positive = [] # List of submitters where ASCAT2, FACETS and Sequenza reported LOH
    posHist = {} # Histogram with the positions (key) and the times a variant is reported in that position (value)
    negative = [] # List of submitters where ASCAT2, FACETS and Sequenza do not report LOH (or less than 2 tools report LOH)
    negHist = {} # Histogram with the positions {key} and times a variant is reported in that position (but for negative submitters)
    pathogenic = [] # List of submitters that are positive for LOH and a pathogenic variant in the gene is found
    patHist = {} # Histogram with the positions {key} and times a variant is reported in that position (for positive submitters with a reported pathogenic variant)
    posData = [] # Two-dimension list with the information of each variant found in LOH submitters
    negData = [] # Two-dimension list with the information of each variant found in no-LOH submitters
    patData = [] # Two-dimension list with the informatino of each variant found in pathogenic submitters

    # Initiallize the data
    clinvarData = readClinVar()
    geneStart = gene[1]
    geneEnd = gene[2]
    for i in range(geneStart, geneEnd+1) :
        posHist[i] = 0
        negHist[i] = 0
        patHist[i] = 0

    print("{} INFO: Getting {} submitters".format(getTime(), cancer))
    # Get the submitter IDs from the cancer repository
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='{cancer}'".format(cancer = cancer))
          cases = q.fetchall()

    print("{} INFO: {} submitters found".format(getTime(), len(cases)))
    # Get the tumors and controls in each submitter
    for c in cases :
        if cases.index(c) % 100 == 0 :
            print("{} INFO: {} submitters analyzed".format(getTime(), cases.index(c)))

        if c[0] not in done : # Don't do the analysis more than once in the same submitter
            # Get the tumor and control samples from the submitter
            with dbcon :
                cur = dbcon.cursor()
                q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
                tumors = q.fetchall()
                q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
                controls = q.fetchall()

            for tm in tumors :
                for cn in controls :
                    pairs += 1
                    print("INFO: {} pairs, {} +, {} -, {} *".format(pairs, len(positive), len(negative), len(pathogenic)))
                    if c[0] not in done :
                        prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
                        # Get and annotate the variants
                        germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
                        cmd = "grep {gene} {vc}".format(gene = genename, vc = germCall)
                        pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                        # Get LOH
                        lohs = 0 # Number of tools that reported LOH (Delecion or copy number normal)
                        loh = []
                        # Get FACETS output file
                        folder = "{wd}/{sub}/{pre}_FACETS".format(wd = wd, sub = c[0], pre = prefix)
                        file = "{fld}/facets_comp_cncf.tsv".format(fld = folder)
                        if os.path.isfile(file) :
                            # Check LOH in the gene, add the output to a list
                            aux = lib.getLOH(file, "facets", gene)
                            loh.append(aux)
                        # Get Sequenza output file
                        folder = "{wd}/{sub}/{pre}_Sequenza".format(wd = wd, sub = c[0], pre = prefix)
                        file = "{fld}/{case}_segments.txt".format(fld = folder, case = c[0])
                        if os.path.isfile(file) :
                            # Check LOH in the gene, add th output to a list
                            aux = lib.getLOH(file, "sequenza", gene)
                            loh.append(aux)
                        # Get PURPLE output file
                        folder = "{wd}/{sub}/{pre}_PURPLE".format(wd = wd, sub = c[0], pre = prefix)
                        file = "{fld}/TUMOR.purple.cnv.somatic.tsv".format(fld = folder)
                        if os.path.isfile(file) :
                            aux = lib.getLOH(file, "purple", gene)
                            loh.append(aux)
                        # Get ASCAT2 output file
                        folder = "{wd}/{sub}/ASCAT2".format(wd = wd, sub = c[0])
                        if os.path.isdir(folder) :
                            # Check LOH in the gene, add the output to a list
                            aux = asc.checkAscat(folder, gene)
                            loh.append(aux)
                        if len(loh) >= 2 :
                            done.append(c[0])
                        # Count the number of LOH found in the patient
                        lohs = loh.count("L") + loh.count("D") # Copy number neutral +`copy number lose`
                        std, err = pr.communicate() # Get the output from grep. This grep was launched before the LOH calling
                        rawVars = std.decode().strip().split("\n")
                        annoVars = []
                        isPathogenic = False
                        for r in rawVars :
                            if r != "" :
                                tmp = annotateClinVar(r, clinvarData)
                                if tmp != "" :
                                    tmp["submitter"] = c[0]
                                    tmp["tumor"] = tm[0]
                                    tmp["control"] = cn[0]
                                    tmp["lohCount"] = str(lohs)
                                    annoVars.append(tmp)
                                    if tmp["type"] in cte.var_positive or tmp["exonicType"] in cte.var_positive :
                                        isPathogenic = True
                                    # # TODO: Tenir en compte que les variants poden estar relacionades amb cancer pero tenir significat com a benign o unknown_significance
                                    if tmp["disease"].find("cancer") > 0 :
                                        isPathogenic = True

                        if lohs >= 2 :
                            positive.append(c[0])
                            # # TODO: Comprovar que les variants stopgain/frameshift es consideren patogeniques tambe
                            if isPathogenic :
                                pathogenic.append(c[0])
                                patData = patData + annoVars
                                for a in annoVars :
                                    aux = int(a["start"])
                                    if aux in patHist :
                                        patHist[aux] += 1
                                    else :
                                        patHist[aux] = 1
                            else :
                                posData = posData + annoVars
                                for a in annoVars :
                                    aux = int(a["start"])
                                    if aux in posHist.keys() :
                                        posHist[aux] += 1
                                    else :
                                        posHist[aux] = 1
                        else :
                            negative.append(c[0])
                            negData += annoVars
                            for a in annoVars :
                                aux = int(a["start"])
                                if aux in negHist.keys() :
                                    negHist[aux] += 1
                                else :
                                    negHist[aux] = 1

    print("{} submitters with enough LOH information".format(len(done)))
    print("{} submitters considered LOH positive".format(len(positive)))
    print("{} submitters had were positive and had a pathogenic variant".format(len(pathogenic)))
    print("{} submitters do not have LOH".format(len(negative)))

    # Save variant position counts in a tsv file
    saveHistogram(posHist, "positiveHistogram.tsv")
    saveHistogram(patHist, "pathogenicHistogram.tsv")
    saveHistogram(negHist, "negativeHistogram.tsv")

    # Save variant data in a tsv file
    saveVariants(posData, "posVariants.tsv")
    saveVariants(patData, "patVariants.tsv")
    saveVariants(negData, "negVariants.tsv")

    return posData, patData, negData

def groupVariants(pos, pat, neg, filename) :
    """Count the number of variants in each group"""
    groups = {}
    addinfo = {}

    for p in pos :
        key = "{chr};{sta};{end};{ref};{alt};{func};{exonic}".format(
            chr = p["chrom"], sta = p["start"], end = p["end"], ref = p["ref"], alt = p["alt"], func = p["type"], exonic = p["exonicType"])
        if key in groups.keys() :
            groups[key]["Positive"] += 1
        else :
            groups[key] = {"Positive" : 1, "Pathogenic" : 0, "Negative" : 0}
            addinfo[key] = {"significance" : p["significance"], "disease" : p["disease"]}

    for p in pat :
        key = "{chr};{sta};{end};{ref};{alt};{func};{exonic}".format(
            chr = p["chrom"], sta = p["start"], end = p["end"], ref = p["ref"], alt = p["alt"], func = p["type"], exonic = p["exonicType"])
        if key in groups.keys() :
            groups[key]["Pathogenic"] += 1
        else :
            groups[key] = {"Positive" : 0, "Pathogenic" : 1, "Negative" : 0}
            addinfo[key] = {"significance" : p["significance"], "disease" : p["disease"]}

    for p in neg :
        key = "{chr};{sta};{end};{ref};{alt};{func};{exonic}".format(
            chr = p["chrom"], sta = p["start"], end = p["end"], ref = p["ref"], alt = p["alt"], func = p["type"], exonic = p["exonicType"])
        if key in groups.keys() :
            groups[key]["Negative"] += 1
        else :
            groups[key] = {"Positive" : 0, "Pathogenic" : 0, "Negative" : 1}
            addinfo[key] = {"significance" : p["significance"], "disease" : p["disease"]}

    with open(filename, "w") as fi :
        fi.write("Chr\tStart\tEnd\tRef\tAlt\tType\tExonicType\tClinVarDisease\tClinVarSignf\tInLOHPositive\tInLOHPathogenic\tInNegative\n")
        for k, v in groups.items() :
            fi.write(k.replace(";", "\t"))
            if k in addinfo :
                fi.write("\t{sig}\t{disease}".format(sig = addinfo[k]["significance"], disease = addinfo[k]["disease"]))
            else :
                fi.write("\tNA\tNA")
            fi.write("\t{}\t".format(v["Positive"]))
            fi.write("\t{}\t".format(v["Pathogenic"]))
            fi.write("{}\n".format(v["Negative"]))

    print("{} INFO: Output stored as {}".format(getTime(), filename))

    return groups

def variantClassifier(vars) :
    p = 0
    n = 0
    u = 0
    signs = {}
    line = ""
    data = {}
    for v in vars :
        # Classify ClinVar significance
        if v["significance"] in signs.keys() :
            signs[v["significance"]] += 1
        else :
            signs[v["significance"]] = 1
        # Classify variants by type
        if v["type"] == "splicing" or (v["type"] == "exonic" and v["exonicType"] in cte.var_positive) :
            p += 1
        elif v["type"] == "exonic" :
            if v["exonicType"] in cte.var_neutral :
                u += 1
            else :
                n += 1
        else :
            n += 1
    data["submitter"] = v["submitter"]
    data["positive"] = p
    data["negative"] = n
    data["unknown"] = u
    for k,v in signs.items() :
        line += "\t{} ({})".format(k, v)
    data["significances"] = line

    return data

def groupSubmitters(variants) :
    current = ""
    tmpVars = []
    allData = []
    for v in variants :
        if current != v["submitter"] :
            current = v["submitter"]
            if len(tmpVars) > 0 :
                allData.append(variantClassifier(tmpVars))
                del(tmpVars)
                tmpVars = [v]
        else :
            tmpVars.append(v)
    return allData

def readFile(path) :
    data = []
    header = True
    with open(path, "r") as fi :
        for l in fi :
            if header :
                header = False
            else :
                aux = l.strip().split("\t")
                dc = {"chrom" : aux[0], "start" : aux[1], "end" : aux[2], "ref" : aux[3], "alt" : aux[4], "type" : aux[5], "exonicType" : aux[6], "lohCount" : aux[7],
                "significance" : aux[8], "disease" : aux[9], "submitter" : aux[10], "tumor" : aux[11], "control" : aux[12]}
                data.append(dc)
    return data


if __name__ == "__main__" :
    # NOTE: To change the analysis parameters, change the constants at the beginning of the file
    # # TODO: If tsv files are created, read them rather than call getData()
    # positive, pathogenic, negative = getData()
    positive = readFile("posVariants.tsv")
    pathogenic = readFile("patVariants.tsv")
    negative = readFile("negVariants.tsv")
    groups = groupVariants(positive, pathogenic, negative, "grouped.vars.tsv")
    # # TODO: Group the variants according to its type. More information in issue #2 on github
    submitters = groupSubmitters(positive)
    tmp = groupSubmitters(pathogenic)
    submitters += tmp
    tmp = groupSubmitters(negative)
    submitters += tmp
    print(len(submitters))
    # # TODO: Run main9.R
