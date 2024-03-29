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
    if d[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] and d[6] == genename:
        anno = []
        chrom = d[0].replace("chr", "")
        # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
        if d[4] == "-" :
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
    if total == 0 :
        print("WARNING: Possible invalid data. No variants found")
        with open(filename, "w") as fi :
            fi.write("position\ttimes\tfreq\n")
            for i in range(minim, maxim+1) :
                if i in data.keys() :
                    fi.write("{}\t{}\t0\n".format(i, data[i]))
                else :
                    fi.write("{}\t0\t0\n".format(i))

        print("{} INFO: Variant histogram data saved as {}".format(getTime(), filename))
    else :
        with open(filename, "w") as fi :
            fi.write("position\ttimes\tfreq\n")
            for i in range(minim, maxim+1) :
                if i in data.keys() :
                    fi.write("{}\t{}\t{}\n".format(i, data[i], round(data[i]/total, 6)))
                else :
                    fi.write("{}\t0\t0\n".format(i))

        print("{} INFO: Variant histogram data saved as {}".format(getTime(), filename))

def saveVariants(data, filename) :
    """Save variants data in a tsv file"""
    with open(filename, "w") as fi :
        fi.write("chrom\tstart\tend\tref\talt\ttype\texonic.type\tlohs\tcln.signf\tcln.disease\tsubmitter\ttumor.id\tcontrol.id\n")
        for d in data :
            fi.write("{cr}\t{st}\t{nd}\t{ref}\t{alt}\t{tp}\t{ex}\t{loh}\t{sgnf}\t{dss}\t{sub}\t{tm}\t{cn}\n".format(
                cr = d["chrom"], st = d["start"], nd = d["end"], ref = d["ref"], alt = d["alt"], tp = d["type"], ex = d["exonicType"], loh = d["lohCount"], sgnf = d["significance"], dss = d["disease"],
                sub = d["submitter"], tm = d["tumor"], cn = d["control"]))

    print("{} INFO: Variants saved as {}".format(getTime(), filename))

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
    positive = {} # List of submitters (and the number of variants in the gene) where ASCAT2, FACETS and Sequenza reported LOH
    posHist = {} # Histogram with the positions (key) and the times a variant is reported in that position (value)
    negative = {} # List of submitters (and the number of variants in the gene) where ASCAT2, FACETS and Sequenza do not report LOH (or less than 2 tools report LOH)
    negHist = {} # Histogram with the positions {key} and times a variant is reported in that position (but for negative submitters)
    pathogenic = {} # List of submitters (and the number of variants in the gene) that are positive for LOH and a pathogenic variant in the gene is found
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
                    if c[0] not in done :
                        prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
                        # Get and annotate the variants
                        germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
                        cmd = "grep -w {gene} {vc}".format(gene = genename, vc = germCall)
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
                            positive[c[0]] = len(annoVars)
                            # # TODO: Comprovar que les variants stopgain/frameshift es consideren patogeniques tambe
                            if isPathogenic :
                                pathogenic[c[0]] = len(annoVars)
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
                            negative[c[0]] = len(annoVars)
                            negData += annoVars
                            for a in annoVars :
                                aux = int(a["start"])
                                if aux in negHist.keys() :
                                    negHist[aux] += 1
                                else :
                                    negHist[aux] = 1

    print("{} INFO: Stats".format(getTime()))
    print("\t{} submitters with enough LOH information".format(len(done)))
    print("\t{} submitters considered LOH positive ({} variants per submitter). {} had no variants in {} gene".format(
        len(positive.keys()), round(sum(positive.values())/len(positive.values()),2), list(positive.values()).count(0), genename))
    print("\t{} submitters considered LOH positive and had a pathogenic variant".format(len(pathogenic.keys())))
    print("\t{} submitters do not have LOH ({} variants per submitter). {} had no variants in {} gene".format(
        len(negative.keys()), round(sum(negative.values())/len(negative.values()),2), list(negative.values()).count(0), genename))

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
        fi.write("Chr\tStart\tEnd\tRef\tAlt\tType\tExonicType\tClinVarDisease\tClinVarSignf\tInLOHPositive\tInLOHPathogenic\tInLOHNegative\n")
        for k, v in groups.items() :
            fi.write(k.replace(";", "\t"))
            if k in addinfo :
                fi.write("\t{disease}\t{sig}".format(sig = addinfo[k]["significance"], disease = addinfo[k]["disease"]))
            else :
                fi.write("\tNA\tNA")
            fi.write("\t{}\t".format(v["Positive"]))
            fi.write("{}\t".format(v["Pathogenic"]))
            fi.write("{}\n".format(v["Negative"]))

    print("{} INFO: Variants grouped and counted by groups. Data stored as {}".format(getTime(), filename))

    return groups

def variantClassifier(vars) :
    """Classify the variant list passed as parameter according to libconstants classification. Additionally, group the Clinvar tags. Return the data in a dictionary"""
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
    """Group the variants by submitter. Classify these variants as positive, negative or unknown. Group ClinVar tags. Return all in a dictionary"""
    tmpVars = []
    allData = []
    if len(variants) > 0 :
        current = variants[0]["submitter"]
        for v in variants :
            if current != v["submitter"] :
                if len(tmpVars) > 0 :
                    allData.append(variantClassifier(tmpVars))
                    del(tmpVars)
                    tmpVars = [v]
                current = v["submitter"]
            else :
                tmpVars.append(v)

        allData.append(variantClassifier(tmpVars))

    return allData

def readFile(path) :
    """Read *Variants.tsv file. Convert the content to the same format as created in getData"""
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

def readCandidateFiles(filename) :
    """Read *candVars.tsv file. Convert the content in a list of dictionaries. Each element in the list is a line (variant) in the file"""
    cand = []
    try :
        with open(filename, "r") as fi :
            header = []
            for l in fi :
                aux = l.strip().split("\t")
                if len(header) == 0 :
                    header = aux
                else :
                    temp = {}
                    it = 0
                    for h in header :
                        temp[h] = aux[it]
                        it += 1
                    cand.append(temp)
    except FileNotFoundError :
        print("WARNING: {} not found".format(filename))

    return cand

def printCandidateVariants(lst) :
    """Print a list of variants (dict) in a fancy format"""
    if len(lst) > 0 :
        for l in lst :
            print("\t{chr}:{st}-{nd}\t{rf}>{lt}_{tp}({ex})".format(chr = l["Chr"], st = l["Start"], nd = l["End"], rf = l["Ref"], lt = l["Alt"], tp = l["Type"], ex = l["ExonicType"]))


def main() :
    """MAIN FUNCTION: Get LOH and variants from a gene, defined as a constant.
    Count how many times a position is mutated in that gene
    Classify the samples (submitters) as
        * Positive: More than one tool reported LOH in that gene
        * Pathogenic: A positive sample that has a pathogenic variant
        * Negative: Only 1 tool (or none) reported LOH in that gene
    Group the variants found, counting how many times the variant is reported in each group
    Count the variants by submitter, classifying the variants as positive, negative or unknown"""
    # NOTE: To change the analysis parameters, change the constants at the beginning of the file
    print("{} INFO: Getting the variants in {} gene written in each {} file".format(getTime(), genename, varCallSuffix))
    # Check if the folder structure has been created before. Ask to remove the data before running the program
    mainFolder = "VUS_project"
    if varCallSuffix == "platypusGerm/platypus.hg38_multianno.txt" :
        subfolder = "{}/Platypus".format(mainFolder)
        geneFolder = "{}/{}".format(subfolder, genename)
    elif varCallSuffix == "strelkaGerm/results/variants/strelka.hg38_multianno.txt" :
        subfolder = "{}/Strelka2".format(mainFolder)
        geneFolder = "{}/{}".format(subfolder, genename)

    if os.path.isdir(geneFolder) :
        opt = input("INPUT: Remove previously calculated data? (Y/N) ")
        if opt == "N" or opt == "n":
            sys.exit()

    # Run the analysis
    # Read the data, either creating the tsv files or reading the previously created files
    positive, pathogenic, negative = getData()
    # positive = readFile("posVariants.tsv")
    # pathogenic = readFile("patVariants.tsv")
    # negative = readFile("negVariants.tsv")

    # Group the variants and count how many times each are present in the different cases
    groups = groupVariants(positive, pathogenic, negative, "grouped.vars.tsv")

    # Group the variants in each submitter according to its type
    submitters = groupSubmitters(positive)
    # Save the data in the last function in a tsv file
    with open("posVarsGrouped.tsv", "w") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in submitters :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))

    tmp = groupSubmitters(pathogenic)
    with open("patVarsGrouped.tsv", "a") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in tmp :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))
    submitters += tmp

    tmp = groupSubmitters(negative)
    with open("negVarsGrouped.tsv", "a") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in tmp :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))
    submitters += tmp

    print("{} INFO: Variant classification grouped by submitter stored as posVarsGrouped.tsv, patVarsGrouped.tsv, and negVarsGrouped.tsv".format(getTime()))

    # Move all the data to a folder
    if not os.path.isdir(mainFolder) :
        os.mkdir(mainFolder)
    if not os.path.isdir(subfolder) :
        os.mkdir(subfolder)
    if not os.path.isdir(geneFolder) :
        os.mkdir(geneFolder)
    os.rename("positiveHistogram.tsv", "{}/positiveHistogram.tsv".format(geneFolder))
    os.rename("pathogenicHistogram.tsv", "{}/pathogenicHistogram.tsv".format(geneFolder))
    os.rename("negativeHistogram.tsv", "{}/negativeHistogram.tsv".format(geneFolder))
    os.rename("posVariants.tsv", "{}/posVariants.tsv".format(geneFolder))
    os.rename("patVariants.tsv", "{}/patVariants.tsv".format(geneFolder))
    os.rename("negVariants.tsv", "{}/negVariants.tsv".format(geneFolder))
    os.rename("grouped.vars.tsv", "{}/grouped.vars.tsv".format(geneFolder))
    os.rename("posVarsGrouped.tsv", "{}/posVarsGrouped.tsv".format(geneFolder))
    os.rename("patVarsGrouped.tsv", "{}/patVarsGrouped.tsv".format(geneFolder))
    os.rename("negVarsGrouped.tsv", "{}/negVarsGrouped.tsv".format(geneFolder))
    print("{temp} INFO: Data moved to {fold}. Run main9.R in {fold} to get further statistics".format(temp = getTime(), fold = geneFolder))

def fastMain() :
    """MAIN FUNCTION: Does the same as main function, but getting the data from a previous main execution. So this function reads posVariants.tsv, patVariants.tsv
    and negVariants.tsv instead of calling getData function"""
    print("{} INFO: Getting the variants and LOH from posVariants, patVariants and negVariants tsv files".format(getTime()))
    positive = readFile("posVariants.tsv")
    pathogenic = readFile("patVariants.tsv")
    negative = readFile("negVariants.tsv")

    # Group the variants and count how many times each are present in the different cases
    groups = groupVariants(positive, pathogenic, negative, "grouped.vars.tsv")

    # Group the variants in each submitter according to its type
    submitters = groupSubmitters(positive)
    # Save the data in the last function in a tsv file
    with open("posVarsGrouped.tsv", "w") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in submitters :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))

    tmp = groupSubmitters(pathogenic)
    with open("patVarsGrouped.tsv", "a") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in tmp :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))
    submitters += tmp

    tmp = groupSubmitters(negative)
    with open("negVarsGrouped.tsv", "a") as fi :
        fi.write("submitter\tvar_positive\tvar_negative\tvar_unknown\tsignificance\n")
        for s in tmp :
            fi.write("{sub}\t{pos}\t{neg}\t{unk}\t{sig}\n".format(sub = s["submitter"], pos = s["positive"], neg = s["negative"], unk = s["unknown"], sig = s["significances"]))
    submitters += tmp

    print("{} INFO: Variant classification grouped by submitter stored as posVarsGrouped.tsv, patVarsGrouped.tsv, and negVarsGrouped.tsv".format(getTime()))

def filterCandidates() :
    """MAIN FUNCTION: Looks for coincidences in Platypus and Strelka2 candidate variants.
    Reports not coincident variants to watch manually what happened
    The coincident variants are searched in cancer repositories where LOH is not expected
    **This function is called after running main9.R in gene folder**
    """
    pl_cand = [] # Candidate variants from Platypus variant caller
    pl_all = [] # All variants from Platypus variant caller
    st_cand = [] # Candidate variants from Strelka2 variant caller
    st_all = [] # All variants from Strelka2 variant caller
    same = []
    dif = []
    abs_path = os.path.abspath(os.getcwd()).split("/")
    gene = abs_path[-1]
    vc = abs_path[-2]
    print("{} INFO: Reading variants in {}".format(getTime(), gene))
    if vc == "Platypus" :
        pl_cand = readCandidateFiles("candVars.tsv")
        st_cand = readCandidateFiles("../../Strelka2/{}/candVars.tsv".format(gene))
        pl_all = readCandidateFiles("allCandVars.tsv")
        st_all = readCandidateFiles("../../Strelka2/{}/allCandVars.tsv".format(gene))
    elif vc == "Strelka2" :
        st_cand = readCandidateFiles("candVars.tsv")
        pl_cand = readCandidateFiles("../../Platypus/{}/candVars.tsv".format(gene))
        st_all = readCandidateFiles("allCandVars.tsv")
        pl_all = readCandidateFiles("../../Platypus/{}/allCandVars.tsv".format(gene))
    else :
        print("ERROR: Variant caller not found")
        sys.exit(1)
    print("{} INFO: {} variants found in Platypus".format(getTime(), len(pl_cand)))
    print("{} INFO: {} variants found in Strelka2".format(getTime(), len(st_cand)))
    for p in pl_cand :
        found = False
        for s in st_cand :
            if p["Start"] == s["Start"] :
                found = True
                break
        if found :
            same.append(p)
        else :
            dif.append(p)
    print("{} INFO: {} variants found in Platypus, but not in Strelka2".format(getTime(), len(dif)))
    printCandidateVariants(dif)
    dif2 = []
    for s in st_cand :
        found = False
        for a in same :
            if a["Start"] == s["Start"] :
                found = True
                break
        if not found :
            dif2.append(s)
    print("{} INFO: {} variants found in Strelka2, but not in Platypus".format(getTime(), len(dif2)))
    printCandidateVariants(dif2)
    print("{} INFO: Searching {} common variants in HNSC. This may take a while".format(getTime(), len(same)))
    print("Gene\tPosition\tType\tReported by\tIn Positive\tIn Pathogenic\tIn Negative\tIn HNSC\tIn Clinvar")
    onlyp = len(pl_cand)
    onlys = len(st_cand)
    both = len(same)
    #Search common variants in HNSC
    for s in same :
        strpth = s["ClinVarSignf"]
        if s["ClinVarSignf"] != "NA" :
            if s["ClinVarSignf"].startswith("Conflicting") :
                strpth = "Conf.Interpret."
            elif s["ClinVarSignf"].startswith("Uncertain") :
                strpth = "VUS"

        if s["Type"] == "exonic" :
            typ = s["ExonicType"]
        else :
            typ = s["Type"]

        cmd = "zgrep -w {coord} /g/strcombio/fsupek_cancer1/TCGA_bam/HNSC/*/*/strelkaGerm/results/variants/variants.vcf.gz | grep -c {chr}".format(coord = s["Start"], chr = s["Chr"])
        pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = pr.communicate()
        print("{gene}\t{pos}\t{type}\tBoth\t{pp}\t{pa}\t{ne}\t{HNSC}\t{pth}".format(
            gene = gene, pos = s["Start"], type = typ, HNSC = out.decode().strip(), pth = strpth, pp = s["InLOHPositive"], pa = s["InLOHPathogenic"], ne = s["InLOHNegative"]))

    # Search only-Platypus variants in HNSC
    for s in dif :
        strpth = s["ClinVarSignf"]
        if s["ClinVarSignf"] != "NA" :
            if s["ClinVarSignf"].startswith("Conflicting") :
                strpth = "Conf.Interpret."
            elif s["ClinVarSignf"].startswith("Uncertain") :
                strpth = "VUS"

        if s["Type"] == "exonic" :
            typ = s["ExonicType"]
        else :
            typ = s["Type"]

        cmd = "zgrep -w {coord} /g/strcombio/fsupek_cancer1/TCGA_bam/HNSC/*/*/strelkaGerm/results/variants/variants.vcf.gz | grep -c {chr}".format(coord = s["Start"], chr = s["Chr"])
        pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = pr.communicate()
        print("{gene}\t{pos}\t{type}\tPlatypus\t{pp}\t{pa}\t{ne}\t{HNSC}\t{pth}".format(
            gene = gene, pos = s["Start"], type = typ, HNSC = out.decode().strip(), pth = strpth, pp = s["InLOHPositive"], pa = s["InLOHPathogenic"], ne = s["InLOHNegative"]))

    # Search only-Strelka2 variants in HNSC
    for s in dif2 :
        strpth = s["ClinVarSignf"]
        if s["ClinVarSignf"] != "NA" :
            if s["ClinVarSignf"].startswith("Conflicting") :
                strpth = "Conf.Interpret."
            elif s["ClinVarSignf"].startswith("Uncertain") :
                strpth = "VUS"

        if s["Type"] == "exonic" :
            typ = s["ExonicType"]
        else :
            typ = s["Type"]

        cmd = "zgrep -w {coord} /g/strcombio/fsupek_cancer1/TCGA_bam/HNSC/*/*/strelkaGerm/results/variants/variants.vcf.gz | grep -c {chr}".format(coord = s["Start"], chr = s["Chr"])
        pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = pr.communicate()
        print("{gene}\t{pos}\t{type}\tStrelka2\t{pp}\t{pa}\t{ne}\t{HNSC}\t{pth}".format(
            gene = gene, pos = s["Start"], type = typ, HNSC = out.decode().strip(), pth = strpth, pp = s["InLOHPositive"], pa = s["InLOHPathogenic"], ne = s["InLOHNegative"]))

if __name__ == "__main__" :
    # Execute main function, depending on the files available
    if os.path.isfile("candVars.tsv") and os.path.isfile("allCandVars.tsv") :
        filterCandidates()
    elif os.path.isfile("posVariants.tsv") and os.path.isfile("patVariants.tsv") and os.path.isfile("negVariants.tsv") :
        fastMain()
    else :
        main()
