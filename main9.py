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
    anno = []
    d = data.split("\t")
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

    return anno


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
                                    aux = a["start"]
                                    patHist[aux] += 1
                            else :
                                posData = posData + annoVars
                                for a in annoVars :
                                    aux = a["start"]
                                    posHist[aux] += 1
                        else :
                            negative.append(c[0])
                            negData += annoVars
                            for a in annoVars :
                                aux = a["start"]
                                negHist += 1

    print("{} submitters with enough LOH information".format(len(done)))
    print("{} submitters considered LOH positive".format(len(positive)))
    print("{} submitters had were positive and had a pathogenic variant".format(len(pathogenic)))
    print("{} submitters do not have LOH".format(len(negative)))
    print(posHist)
    
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
                fi.write("\t")
                fi.write(addinfo[k])
            else :
                fi.write("\tNA\tNA")
            fi.write("\t{}\t".format(v["Positive"]))
            fi.write("\t{}\t".format(v["Pathogenic"]))
            fi.write("{}\n".format(v["Negative"]))

    print("{} INFO: Output stored as {}".format(getTime(), filename))

    return groups


if __name__ == "__main__" :
    # NOTE: To change the analysis parameters, change the constants at the beginning of the file
    positive, pathogenic, negative = getData()
    groups = groupVariants(positive, pathogenic, negative, "grouped.vars.tsv")



#
# def getData() :
#     """Get LOH and variants information for each submitter in the gene of interest
#
#     Classify the submitters as
#         * Positive: When 2 or more LOH tools report LOH (CNN-LOH or CNL-LOH) in the gene of interest
#         * Negative: When less than 2 LOH tools report LOH in the gene of interest
#
#     Variants are classified in each list (posData or negData) according to LOH found
#     """
#     # Variables
#     pairs = 0 # Number of submitters with pair tumor-control
#     done = [] # List of submitters with ASCAT2, FACETS and Sequenza done
#     positive = [] # List of submitters where ASCAT2, FACETS and Sequenza reported LOH
#     posHist = {} # Histogram with the positions (key) and the times a variant is reported in that position (value)
#     negative = [] # List of submitters where ASCAT2, FACETS and Sequenza do not report LOH (or less than 2 tools report LOH)
#     negHist = {} # Histogram with the positions {key} and times a variant is reported in that position (but for negative submitters)
#     posData = [] # Two-dimension list with the information of each variant found in LOH submitters
#     negData = [] # Two-dimension list with the information of each variant found in no-LOH submitters
#
#     print("{} INFO: Getting {} submitters".format(getTime(), cancer))
#     # Get the submitter IDs from the cancer repository
#     with dbcon :
#           cur = dbcon.cursor()
#           q = cur.execute("SELECT submitter FROM patient WHERE cancer='{cancer}'".format(cancer = cancer))
#           cases = q.fetchall()
#
#     print("{} INFO: {} submitters found".format(getTime(), len(cases)))
#     # Get the tumors and controls in each submitter
#     for c in cases :
#         if cases.index(c) % 100 == 0 :
#             print("{} INFO: {} submitters analyzed".format(getTime(), cases.index(c)))
#
#         if c[0] not in done : # Don't do the analysis more than once in the same submitter
#             # Get the tumor and control samples from the submitter
#             with dbcon :
#                 cur = dbcon.cursor()
#                 q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
#                 tumors = q.fetchall()
#                 q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
#                 controls = q.fetchall()
#
#             for tm in tumors :
#                 for cn in controls :
#                     pairs += 1
#                     loh = [] # Store the LOH reported by the tools
#                     prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
#                     # Get FACETS output file
#                     folder = "{wd}/{sub}/{pre}_FACETS".format(wd = wd, sub = c[0], pre = prefix)
#                     file = "{fld}/facets_comp_cncf.tsv".format(fld = folder)
#                     if os.path.isfile(file) :
#                         # Check LOH in the gene, add the output to a list
#                         aux = lib.getLOH(file, "facets", gene)
#                         loh.append(aux)
#                     # Get Sequenza output file
#                     folder = "{wd}/{sub}/{pre}_Sequenza".format(wd = wd, sub = c[0], pre = prefix)
#                     file = "{fld}/{case}_segments.txt".format(fld = folder, case = c[0])
#                     if os.path.isfile(file) :
#                         # Check LOH in the gene, add th output to a list
#                         aux = lib.getLOH(file, "sequenza", gene)
#                         loh.append(aux)
#                     # Get PURPLE output file
#                     folder = "{wd}/{sub}/{pre}_PURPLE".format(wd = wd, sub = c[0], pre = prefix)
#                     file = "{fld}/TUMOR.purple.cnv.somatic.tsv".format(fld = folder)
#                     if os.path.isfile(file) :
#                         aux = lib.getLOH(file, "purple", gene)
#                         loh.append(aux)
#                     # Get ASCAT2 output file
#                     folder = "{wd}/{sub}/ASCAT2".format(wd = wd, sub = c[0])
#                     if os.path.isdir(folder) :
#                         # Check LOH in the gene, add the output to a list
#                         aux = asc.checkAscat(folder, gene)
#                         loh.append(aux)
#                     if len(loh) >= 2 :
#                         done.append(c[0])
#                     # Count the number of LOH found in the patient
#                     lohs = loh.count("L") + loh.count("D") # Copy number neutral +`copy number lose`
#                     if lohs >= 2 :
#                         positive.append(c[0])
#                         # Get the gene variants
#                         germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
#                         cmd = "grep {gene} {vc}".format(vc = germCall, gene = genename)
#                         pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#                         std, err = pr.communicate()
#                         out = std.decode().strip().split("\n")
#                         for o in out :
#                             aux = o.split("\t")
#                             if len(aux) > 1 :
#                                 # Create a histogram that counts the frequency of each position
#                                 pos = int(aux[1])
#                                 if pos in posHist.keys()  :
#                                     posHist[pos] += 1
#                                 else :
#                                     posHist[pos] = 1
#                                 # In case we prefer to collect the full variant...
#                                 # key = "{}-{}-{}".format(pos, ref, alt)
#                                 # if key in variants.keys() :
#                                 #     variants[key] += 1
#                                 # else :
#                                 #     variants[key] = 1
#                                 # Store the variant information in a variable
#                                 # Do not add intergenic, downstream, upstream variants
#                                 if aux[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] :
#                                     # Column order: chrom, start, end, ref, alt, gene_position, exonic_type, submitter, uuid_tumor, uuid_control, lohs_reported
#                                     posData.append([aux[0], aux[1], aux[2], aux[3], aux[4], aux[5], aux[8], c[0], tm[0], cn[0], str(lohs)])
#
#                     else :
#                         negative.append(c[0])
#                         germCall = "{wd}/{sub}/{uuid}/{suffix}".format(wd = wd, sub = c[0], uuid = cn[0], suffix = varCallSuffix)
#                         cmd = "grep {gene} {vc}".format(vc = germCall, gene = genename)
#                         pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#                         std, err = pr.communicate()
#                         out = std.decode().strip().split("\n")
#                         for o in out :
#                             aux = o.split("\t")
#                             if len(aux) > 1 :
#                                 # Create a histogram that counts the frequency of each position
#                                 pos = int(aux[1])
#                                 if pos in negHist.keys()  :
#                                     negHist[pos] += 1
#                                 else :
#                                     negHist[pos] = 1
#                                 # Do not add intergenic, downstream, upstream variants
#                                 if aux[5] in ["intronic", "exonic", "splicing", "UTR3", "UTR5"] :
#                                     # Column order: chrom, start, end, ref, alt, gene_position, exonic_type, submitter, uuid_tumor, uuid_control, lohs_reported
#                                     negData.append([aux[0], aux[1], aux[2], aux[3], aux[4], aux[5], aux[8], c[0], tm[0], cn[0], str(lohs)])
#
#
#     print("{} INFO: Pairs tumor-control: {}".format(getTime(), pairs))
#     print("{} INFO: Analysis done in {} pairs".format(getTime(), len(done)))
#     print("{} INFO: {}  had 2 or more LOH reported in {} gene".format(getTime(), len(positive), genename))
#
#     # Print the data in histogram format to do a plot in R that searches for clusters
#     # Two possible options to plot the variants
#     # 1) From the lowest mutation coordinate
#     # minim = min(posHist.keys())
#     # maxim = max(posHist.keys())
#     # 2) From the gene start
#     minim = gene[1]
#     maxim = gene[2]
#     with open("positionHistogram.tsv", "w") as fi :
#         fi.write("position\ttimes\n")
#         for i in range(minim, maxim+1) :
#             if i in posHist.keys() :
#                 fi.write("{}\t{}\n".format(i, posHist[i]))
#             else :
#                 fi.write("{}\t0\n".format(i))
#
#     with open("negativeHistogram.tsv", "w") as fi :
#         fi.write("position\ttimes\n")
#         for i in range(minim, maxim+1) :
#             if i in negHist.keys() :
#                 fi.write("{}\t{}\n".format(i, negHist[i]))
#             else :
#                 fi.write("{}\t0\n".format(i))
#
#     print("{} INFO: Variant-position histograms stored as positionHistogram.tsv and negativeHistogram.tsv".format(getTime()))
#
#     # Write the variants data in a tab-separated file
#     with open("posVariants.tsv", "w") as fi :
#         for l in posData :
#             fi.write("\t".join(l))
#             fi.write("\n")
#
#     with open("negVariants.tsv", "w") as fi :
#         for l in negData :
#             fi.write("\t".join(l))
#             fi.write("\n")
#
#     print("{} INFO: Variant information stored as posVariants.tsv and negVariants.tsv".format(getTime()))
#     print("{} INFO: Creating the plots".format(getTime()))
#     cmd = "Rscript main9.R {}".format(genename)
#     pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     std, err = pr.communicate()
#     print(std.decode())
#     if pr.returncode != 0 :
#         print("{} WARNING: Error found while running R ({}). Description\n".format(getTime(), cmd, err.decode()))
#
#     # Return variants data for further analysis
#     return posData, negData
#
#
#
# def filterVariants(data, filename, cnt = None) :
#     """Annotate the variants passed as parameter (data) with ClinVar data (cnt variable. If data is not available yet, it calls the function to get ClinVar information)
#     Then reads the variant information obtained and classifies the submitters according to the variant pathogenicity:
#         * Pathogenic when a pathogenic variant (splicing, exonic-indels, exonic-stopgain) is found. Or when the variant is annotated in ClinVar as related with cancer
#         * Negative when no pathogenic variant is found, neither the variants are annotated as related with cancer
#     """
#     posSubmitters = [] # Submitters where a pathogenic variant is found
#     tmplist = []
#     nopos = [] # Variants that are in submitters with a pathogenic variant (like frameshift indels, stopgains...)
#     ispos = []
#
#     # Get ClinVar data
#     if cnt == None :
#         cnt = readClinVar()
#
#     for d in data :
#         if len(d) > 8 :
#             # Check if the variant type is considered positive
#             if d[6] in cte.var_positive or d[5] in cte.var_positive:
#                 posSubmitters.append(d[8])
#
#             # As ANNOVAR changes the position in insertions/deletions, we substract 1 to the start position
#             if d[4] == "-" :
#                 pos = int(d[1]) - 1
#             else :
#                 pos = int(d[1])
#
#             search = "{}-{}".format(d[0].replace("chr", ""), pos)
#             supData = {"db" : {}, "disease" : "NA", "significance" : "NA", "revStatus" : "NA"}
#             if search in cnt.keys() :
#                 supData = getClinVar(cnt[search], d)
#
#             # Check if the variant is associated with cancer
#             if supData["disease"].find("cancer") > 0 and d[8] not in posSubmitters :
#                 posSubmitters.append(d[8])
#
#             # Check if the submitter is considered a positive case
#             d.append(supData["disease"])
#             d.append(supData["significance"])
#             tmplist.append(d)
#
#     for v in tmplist :
#         if v[8] in posSubmitters :
#             ispos.append(v)
#         else :
#             nopos.append(v)
#
#     print("{} INFO: {} submitters removed".format(getTime(), len(posSubmitters)))
#     print("{} INFO: {} variants removed. Data saved as {}_pathogenic.tsv".format(getTime(), len(ispos), filename))
#     with open("{}_pathogenic.tsv".format(filename), "w") as fi :
#         for v in ispos :
#             fi.write("\t".join(v))
#             fi.write("\n")
#     print("{} INFO: {} variants preserved. Data saved as {}_negative.tsv".format(getTime(), len(nopos), filename))
#     with open("{}_negative.tsv".format(filename), "w") as fi :
#         for v in nopos :
#             fi.write("\t".join(v))
#             fi.write("\n")
#
#     return ispos, nopos
#
#
# if __name__ == "__main__" :
#     print("{} INFO: Getting the variants in {} gene written in each {} file".format(getTime(), genename, varCallSuffix))
#     pos_variants, neg_variants = getData()
#     clinvar = readClinVar()
#     # TODO: MODIFICAR AQUESTES LINIES PER AGAFAR LES VARIANTS EN EL NOU FORMAT
#     # with open("posVariants.tsv", "r") as fi :
#     #     pos_variants = fi.read()
#     # with open("negVariants.tsv", "r") as fi :
#     #     neg_variants = fi.read()
#     patho, nega = filterVariants(pos_variants, "posVariants.annotated", clinvar)
#     groupVariants(patho, nega, "posVariants.grouped.tsv")
#     # patho, nega = filterVariants(neg_variants, "negVariants.annotated", clinvar)
#     # groupVariants(patho, nega, "negVariants.grouped.tsv")
#     # aux = []
#     # for d in pos_variants.split("\n") :
#     #     if d != "" :
#     #         aux.append(d.split("\t"))
#     # pos_variants = aux
#     # aux = []
#     # for d in neg_variants.split("\n") :
#     #     if d != "" :
#     #         aux.append(d.split("\t"))
#     # neg_variants = aux
#     # pos_variants = annotateClinVar(pos_variants, clinvar)
#     # groupVariants(pos_variants, neg_variants, "allVariants.grouped.tsv")
