#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import subprocess
import sys

import main1 as lib
import main8 as asc

# Constants
brca1 = ["17", 43044295, 43125483] # Gene of interest to search LOH
cFolder = "fsupek_cancer2"
cancer = "OV" # Cancer repository to search
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")

wd = "/g/strcombio/{cancer_path}/TCGA_bam/{c}".format(c = cancer, cancer_path = cFolder)

# Stats
pairs = 0 # Number of submitters with pair tumor-control
done = [] # List of submitters with ASCAT2, FACETS and Sequenza done
positive = [] # List of submitters where ASCAT2, FACETS and Sequenza reported LOH
variants = {} # Histogram with the positions (key) and the times a variant is reported in that position (value)
negative = [] # List of submitters where ASCAT2, FACETS and Sequenza do not report LOH (or less than 2 tools report LOH)
negVars = {} # Histogram with the positions {key} and times a variant is reported in that position (but for negative submitters)
data = "" # Text table with the information about the variants found
negData = ""

# In the case we want to do it more automatically
gene = brca1
genename = "BRCA1"
varCallSuffix = "platypusGerm/platypus.hg38_multianno.txt"
#varCallSuffix = "strelkaGerm/results/variants/strelka.hg38_multianno.txt"

# Get the submitter IDs from the cancer repository
with dbcon :
      cur = dbcon.cursor()
      q = cur.execute("SELECT submitter FROM patient WHERE cancer='{cancer}'".format(cancer = cancer))
      cases = q.fetchall()

print("INFO: Total submitters in {} cancer: {}".format(cancer, len(cases)))
# Get the tumors and controls in each submitter
for c in cases :
    if cases.index(c) % 100 == 0 :
        print("INFO: {} cases executed".format(cases.index(c)))
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
            # Get ASCAT2 output file
            folder = "{wd}/{sub}/ASCAT2".format(wd = wd, sub = c[0])
            if os.path.isdir(folder) :
                # Check LOH in the gene, add the output to a list
                aux = asc.checkAscat(folder, gene)
                loh.append(aux)
            if len(loh) == 3 :
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
                        ref = aux[2]
                        alt = aux[3]
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
                        negData += "{chr}\t{st}\t{end}\t{ref}\t{alt}\t{ex}\t{typex}\t{sub}\t{idtm}\t{idcn}\n".format(chr = aux[0], st = aux[1], end = aux[2], ref = aux[3], alt = aux[4], ex = aux[5], typex = aux[8], sub = c[0], idtm = tm[0], idcn = cn[0])


# Check if FACETS/Sequenza/ASCAT2 have reported LOH
print("INFO: Pairs tumor-control: {}".format(pairs))
print("INFO: Analysis done in {} pairs".format(len(done)))
print("INFO: {}  had 2 or more LOH reported in {} gene".format(len(positive), genename))

print("INFO: Positive cases. Variant classification")
sortList = sorted(variants, key = variants.get, reverse = True)
for k in sortList :
    print("{} - {}".format(variants[k], k))

# Print the data in histogram format to do a plot in R that searches for clusters
# Two possible options to plot the variants
# 1) From the lowest mutation coordinate
minim = min(variants.keys())
maxim = max(variants.keys())
# 2) From the gene start
minim = brca1[1]
maxim = brca1[2]
with open("positionHistogram.tsv", "w") as fi :
    fi.write("position\ttimes\n")
    for i in range(minim, maxim+1) :
        if i in variants.keys() :
            fi.write("{}\t{}\n".format(i, variants[i]))
        else :
            fi.write("{}\t0\n".format(i))

print("INFO: Negative cases. Variant classification")
sortList = sorted(negVars, key = negVars.get, reverse = True)
for k in sortList :
    print("{} - {}".format(negVars[k], k))

with open("negativeHistogram.tsv", "w") as fi :
    fi.write("position\ttimes\n")
    for i in range(minim, maxim+1) :
        if i in negVars.keys() :
            fi.write("{}\t{}\n".format(i, negVars[i]))
        else :
            fi.write("{}\t0\n".format(i))

print("INFO: Variant-position histogram stored as positionHistogram.tsv and negativeHistogram.tsv")

with open("posVariants.tsv", "w") as fi :
    fi.write(data)

with open("negVariants.tsv", "w") as fi :
    fi.write(negData)
print("INFO: Variant information stored as variants.tsv")
