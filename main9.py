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

# Constants
brca1 = ["17", 43044295, 43125483] # Gene of interest to search LOH
brca2 = ["13", 32315508, 32400268]
atm = ["11", 108222832, 108369099]
palb2 = ["16", 23603165, 23641310]
cFolder = "fsupek_cancer2"
cancer = "OV" # Cancer repository to search
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/{cancer_path}/TCGA_bam/{c}".format(c = cancer, cancer_path = cFolder)

# In the case we want to do it more automatically
gene = brca1
genename = "BRCA1"
varCallSuffix = "platypusGerm/platypus.hg38_multianno.txt"
#varCallSuffix = "strelkaGerm/results/variants/strelka.hg38_multianno.txt"

def getTime() :
    """Returns current time (hours:minutes:seconds) in a fancy format"""
    return time.strftime("%H:%M:%S", time.gmtime())

def getClinVar(clinvar, variant) :
    """Extract the important information from ClinVar

    Columns extracted (description given in clinvar vcf file)
        * CLNDISDB: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
        * CLNDN: ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
        * CLNSIG: Clinical significance for this single variant
        * CLNREVSTAT: ClinVar review status for the Variation ID
    """

    data = {"db" : "NA", "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
    # Check if the variant is exactly the same as the reported by ClinVar
    tmp = variant.split("\t")
    ref = tmp[3]
    alt = tmp[4]
    aux = tmp[0:5]
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

def getData() :
    # Variables
    # Stats
    pairs = 0 # Number of submitters with pair tumor-control
    done = [] # List of submitters with ASCAT2, FACETS and Sequenza done
    positive = [] # List of submitters where ASCAT2, FACETS and Sequenza reported LOH
    variants = {} # Histogram with the positions (key) and the times a variant is reported in that position (value)
    negative = [] # List of submitters where ASCAT2, FACETS and Sequenza do not report LOH (or less than 2 tools report LOH)
    negVars = {} # Histogram with the positions {key} and times a variant is reported in that position (but for negative submitters)
    data = "" # Text table with the information about the variants found
    negData = ""

    print("{} INFO: Getting {} submitters".format(getTime(), cancer))
    # Get the submitter IDs from the cancer repository
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='{cancer}'".format(cancer = cancer))
          cases = q.fetchall()

    print("{} INFO: Total submitters in {} cancer: {}".format(getTime(), cancer, len(cases)))
    # Get the tumors and controls in each submitter
    for c in cases :
        if cases.index(c) % 100 == 0 :
            print("{} INFO: {} cases executed".format(getTime(), cases.index(c)))

        with dbcon :
            cur = dbcon.cursor()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
            tumors = q.fetchall()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
            controls = q.fetchall()

        for tm in tumors :
            if c[0] in done : # Don't do the analysis more than once in the same submitter
                break
            for cn in controls :
                pairs += 1
                loh = [] # Store the LOH reported by the tools
                prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
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
                if len(loh) >= 3 :
                    done.append(c[0])
                # Count the number of LOH found in the patient
                lohs = loh.count("L") + loh.count("D") # Copy number neutral +`copy number lose`
                if lohs >= 2 :
                    positive.append(c[0])
                    # Get the gene variants
                    germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
                    cmd = "grep {gene} {vc}".format(vc = germCall, gene = genename)
                    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    std, err = pr.communicate()
                    out = std.decode().strip().split("\n")
                    for o in out :
                        aux = o.split("\t")
                        if len(aux) > 1 :
                            # Create a histogram that counts the frequency of each position
                            pos = int(aux[1])
                            if pos in variants.keys()  :
                                variants[pos] += 1
                            else :
                                variants[pos] = 1
                            # In case we prefer to collect the full variant...
                            # key = "{}-{}-{}".format(pos, ref, alt)
                            # if key in variants.keys() :
                            #     variants[key] += 1
                            # else :
                            #     variants[key] = 1
                            # Store the variant information in a variable
                            # Do not add intergenic, downstream, upstream variants
                            if aux[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] :
                                data += "{chr}\t{st}\t{end}\t{ref}\t{alt}\t{ex}\t{typex}\t{sub}\t{idtm}\t{idcn}\n".format(chr = aux[0], st = aux[1], end = aux[2], ref = aux[3], alt = aux[4], ex = aux[5], typex = aux[8], sub = c[0], idtm = tm[0], idcn = cn[0])
                else :
                    negative.append(c[0])
                    germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
                    cmd = "grep {gene} {vc}".format(vc = germCall, gene = genename)
                    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    std, err = pr.communicate()
                    out = std.decode().strip().split("\n")
                    for o in out :
                        aux = o.split("\t")
                        if len(aux) > 1 :
                            # Create a histogram that counts the frequency of each position
                            pos = int(aux[1])
                            ref = aux[2]
                            alt = aux[3]
                            if pos in negVars.keys()  :
                                negVars[pos] += 1
                            else :
                                negVars[pos] = 1
                            # Do not add intergenic, downstream, upstream variants
                            if aux[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] :
                                negData += "{chr}\t{st}\t{end}\t{ref}\t{alt}\t{ex}\t{typex}\t{sub}\t{idtm}\t{idcn}\n".format(chr = aux[0], st = aux[1], end = aux[2], ref = aux[3], alt = aux[4], ex = aux[5], typex = aux[8], sub = c[0], idtm = tm[0], idcn = cn[0])


    # Check if FACETS/Sequenza/ASCAT2 have reported LOH
    print("{} INFO: Pairs tumor-control: {}".format(getTime(), pairs))
    print("{} INFO: Analysis done in {} pairs".format(getTime(), len(done)))
    print("{} INFO: {}  had 2 or more LOH reported in {} gene".format(getTime(), len(positive), genename))

    # Print the data in histogram format to do a plot in R that searches for clusters
    # Two possible options to plot the variants
    # 1) From the lowest mutation coordinate
    minim = min(variants.keys())
    maxim = max(variants.keys())
    # 2) From the gene start
    minim = gene[1]
    maxim = gene[2]
    with open("positionHistogram.tsv", "w") as fi :
        fi.write("position\ttimes\n")
        for i in range(minim, maxim+1) :
            if i in variants.keys() :
                fi.write("{}\t{}\n".format(i, variants[i]))
            else :
                fi.write("{}\t0\n".format(i))

    with open("negativeHistogram.tsv", "w") as fi :
        fi.write("position\ttimes\n")
        for i in range(minim, maxim+1) :
            if i in negVars.keys() :
                fi.write("{}\t{}\n".format(i, negVars[i]))
            else :
                fi.write("{}\t0\n".format(i))

    print("{} INFO: Variant-position histogram stored as positionHistogram.tsv and negativeHistogram.tsv".format(getTime()))

    with open("posVariants.tsv", "w") as fi :
        fi.write(data)

    with open("negVariants.tsv", "w") as fi :
        fi.write(negData)

    print("{} INFO: Variant information stored as variants.tsv".format(getTime()))
    print("{} INFO: Creating the plots".format(getTime()))
    cmd = "Rscript main9.R {}".format(genename)
    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = pr.communicate()
    print(std.decode())
    if pr.returncode != 0 :
        print("WARNING: Error found while running R. Description\n".format(err.decode()))

    return data, negData

def filterVariants(data, filename) :
    """Check the variants and remove the submitters (and their variants) with a pathogenic (positive) variant"""
    posSubmitters = [] # Submitters where a pathogenic variant is found
    nopos = [] # Variants that are not considered positive (like frameshift indels, stopgains...)
    ispos = []
    keyword = "cancer" # If the variant has this keyword. It can be considered pathogenic
    vcf = "clinvar.vcf"
    cnt = {} # Data from clinvar vcf
    # Get ClinVar data
    with open(vcf, "r") as fi :
        for l in fi :
            if not l.startswith("#") :
                aux = l.strip().split("\t")
                idx = "{}-{}".format(aux[0], aux[1])
                if idx not in cnt.keys() : # PATCH!! Get the first element if there are position duplicates
                    cnt[idx] = l

    for d in data.split("\n") :
        v = d.split("\t")
        if len(v) > 8 :
            # Check if the variant type is considered positive
            if v[6] in cte.var_positive :
                posSubmitters.append(v[7])

            # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
            if aux[4] == "-" :
                pos = int(aux[1]) - 1
            else :
                pos = int(aux[1])

            search = "{}-{}".format(aux[0].replace("chr", ""), pos)
            supData = {"db" : {}, "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
            if search in cnt.keys() :
                supData = getClinVar(cnt[search], v)

            # Check if the variant is associated with cancer
            if supData["disease"].find("cancer") > 0 and aux[7] not in posSubmitters :
                posSubmitters.append(aux[7])

            # Check if the submitter is considered a positive case
            aux.append(supData["disease"])
            aux.append(supData["significance"])
            if v[7] in posSubmitters :
                ispos.append(aux)
            else :
                nopos.append(aux)

    print("INFO: {} submitters removed".format(len(posSubmitters)))
    print("INFO: {} variants removed".format(len(ispos)))
    print("INFO: {} variants preserved".format(len(nopos)))



def old_filterVariants(data, filename) :
    # Post-production. Check the positive variants and remove the submitters with a pathogenic variant
    posSubmitters = [] # Submitters with a pathogenic variant found
    posData = []
    txt = ""
    for d in data.split("\n") :
        v = d.split("\t")
        if len(v) > 8 :
            if v[7] not in posSubmitters : #Only check variants that are not in submitters with pathogenic variant
                if v[6] in cte.var_positive :
                    posSubmitters.append(v[7])
                else :
                    posData.append(d)

    print("INFO: Removed {} submitters with a pathogenic variant".format(len(posSubmitters)))
    print("INFO: Removed {} variants".format(len(data.split("\n")) - len(posData)))

    vcf = "clinvar.vcf"
    cnt = {}
    with open(vcf, "r") as fi :
        for l in fi :
            if not l.startswith("#") :
                aux = l.strip().split("\t")
                idx = "{}-{}".format(aux[0], aux[1])
                if idx not in cnt.keys() : # PATCH!! Get the first element if there are position duplicates
                    cnt[idx] = l

    # Check if the variant is reported in ClinVar
    keyword = "cancer" # If the variant has this keyword. It can be considered pathogenic
    for v in posData :
        aux = v.split("\t")
        # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
        if aux[4] == "-" :
            pos = int(aux[1]) - 1
        else :
            pos = int(aux[1])

        search = "{}-{}".format(aux[0].replace("chr", ""), pos)
        supData = {"db" : {}, "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
        if search in cnt.keys() :
            supData = getClinVar(cnt[search], v)

        if supData["disease"].find("cancer") > 0 and aux[7] not in posSubmitters :
            posSubmitters.append(aux[7])

        txt += "{data}\t{dss}\t{sig}\n".format(data = v.strip(), dss = supData["disease"], sig = supData["significance"])

    print("INFO: Removed {} submitters".format(len(posSubmitters)))
    auxList = []
    # Remove the variants that come from a positive submitter (a submitter with a positive/pathogenic variant)
    for v in posData :
        aux = v.split("\t")
        if aux[7] not in posSubmitters :
            auxList.append(v)

    print("INFO: Removed {} variants".format(len(posData) - len(auxList)))
    posData = auxList

    print("{} INFO: ClinVar annotated variants stored as {}".format(getTime(), filename))
    with open(filename, "w") as fi :
        fi.write(txt)



# # Create a new dict with the variants from the submitters that are not considered positive
# posVars = {}
# for k, v in variants.items() :
#     if v[7] not in posSubmitters :
#         posVars

if __name__ == "__main__" :
    # pos_variants, neg_variants = getData()
    with open("posVariants.tsv", "r") as fi :
        pos_variants = fi.read()
    filterVariants(pos_variants, "posVariants.annotated.tsv")
    filterVariants(neg_variants, "negVariants.annotated.tsv")
