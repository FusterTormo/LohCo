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
gene = palb2
genename = "PALB2"
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
                                    data += "{chr}\t{st}\t{end}\t{ref}\t{alt}\t{ex}\t{typex}\t{loh}\t{sub}\t{idtm}\t{idcn}\n".format(
                                        chr = aux[0], st = aux[1], end = aux[2], ref = aux[3], alt = aux[4], ex = aux[5], typex = aux[8], sub = c[0], idtm = tm[0], idcn = cn[0], loh = lohs)
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
                                    negData += "{chr}\t{st}\t{end}\t{ref}\t{alt}\t{ex}\t{typex}\t{loh}\t{sub}\t{idtm}\t{idcn}\n".format(
                                        chr = aux[0], st = aux[1], end = aux[2], ref = aux[3], alt = aux[4], ex = aux[5], typex = aux[8], sub = c[0], idtm = tm[0], idcn = cn[0], loh = lohs)


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
        print("{} WARNING: Error found while running R. Description\n".format(getTime(), err.decode()))

    return data, negData

def filterVariants(data, filename) :
    """Check the variants and remove the submitters (and their variants) with a pathogenic (positive) variant"""
    posSubmitters = [] # Submitters where a pathogenic variant is found
    tmplist = []
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
            if v[6] in cte.var_positive or v[5] in cte.var_positive:
                posSubmitters.append(v[8])

            # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
            if v[4] == "-" :
                pos = int(v[1]) - 1
            else :
                pos = int(v[1])

            search = "{}-{}".format(v[0].replace("chr", ""), pos)
            supData = {"db" : {}, "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
            if search in cnt.keys() :
                supData = getClinVar(cnt[search], d)

            # Check if the variant is associated with cancer
            if supData["disease"].find("cancer") > 0 and v[8] not in posSubmitters :
                posSubmitters.append(v[8])

            # Check if the submitter is considered a positive case
            v.append(supData["disease"])
            v.append(supData["significance"])
            tmplist.append(v)

    for v in tmplist :
        if v[8] in posSubmitters :
            ispos.append(v)
        else :
            nopos.append(v)

    print("{} INFO: {} submitters removed".format(getTime(), len(posSubmitters)))
    print("{} INFO: {} variants removed. Data saved as {}_pathogenic.tsv".format(getTime(), len(ispos), filename))
    with open("{}_pathogenic.tsv".format(filename), "w") as fi :
        for v in ispos :
            fi.write("\t".join(v))
            fi.write("\n")
    print("{} INFO: {} variants preserved. Data saved as {}_negative.tsv".format(getTime(), len(nopos), filename))
    with open("{}_negative.tsv".format(filename), "w") as fi :
        for v in nopos :
            fi.write("\t".join(v))
            fi.write("\n")

    return ispos, nopos

def groupVariants(patho, nega, filename) :
    """Count the number of variants in each group"""
    groups = {}

    for n in nega :
        key = "{chr};{sta};{end};{ref};{alt};{func};{exonic}".format(chr = n[0], sta = n[1], end = n[2], ref = n[3], alt = n[4], func = n[5], exonic = n[6])
        if key in groups.keys() :
            groups[key]["No_pathogenic"] += 1
        else :
            groups[key] = {"No_pathogenic" : 1, "Pathogenic" : 0}

    for p in patho :
        key = "{chr};{sta};{end};{ref};{alt};{func};{exonic}".format(chr = p[0], sta = p[1], end = p[2], ref = p[3], alt = p[4], func = p[5], exonic = p[6])
        if key in groups.keys() :
            groups[key]["Pathogenic"] += 1
        else :
            groups[key] = {"Pathogenic" : 1, "No_pathogenic" : 0}

    with open(filename, "w") as fi :
        fi.write("Chr\tStart\tEnd\tRef\tAlt\tType\tExonicType\tInNegative\tInPathogenic\n")
        for k, v in groups.items() :
            fi.write(k.replace(";", "\t"))
            fi.write("\t{}\t".format(v["No_pathogenic"]))
            fi.write("{}\n".format(v["Pathogenic"]))

    print("{} INFO: Output stored as {}".format(getTime(), filename))

if __name__ == "__main__" :
    print("{} INFO: Getting the variants in {} gene written in each {} file".format(getTime(), genename, varCallSuffix))
    pos_variants, neg_variants = getData()
    # with open("posVariants.tsv", "r") as fi :
    #     pos_variants = fi.read()
    # with open("negVariants.tsv", "r") as fi :
    #     neg_variants = fi.read()
    patho, nega = filterVariants(pos_variants, "posVariants.annotated")
    groupVariants(patho, nega, "posVariants.grouped.tsv")
    patho, nega = filterVariants(neg_variants, "negVariants.annotated")
    groupVariants(patho, nega, "negVariants.grouped.tsv")
    aux = []
    for d in pos_variants.split("\n") :
        if d != "" :
            aux.append(d.split("\t"))
    pos_variants = aux
    aux = []
    for d in neg_variants.split("\n") :
        if d != "" :
            aux.append(d.split("\t"))
    neg_variants = aux
    groupVariants(pos_variants, neg_variants, "allVariants.grouped.tsv")
