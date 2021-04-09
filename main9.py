#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import subprocess
import sys

import main1 as lib
import main8 as asc

# Constants
brca1 = ["17", 43044295, 43170245] # Gene of interest to search LOH
cFolder = "fsupek_cancer2"
cancer = "OV" # Cancer repository to search
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")

wd = "/g/strcombio/{cancer_path}/TCGA_bam/{c}".format(c = cancer, cancer_path = cFolder)

# Stats
pairs = 0
done = []
positive = []

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

print("INFO: Total submitters: {}".format(len(cases)))
# Get the tumors and controls in each submitter
for c in cases :
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
                    



# Check if FACETS/Sequenza/ASCAT2 have reported LOH
print("INFO: Pairs tumor-control: {}".format(pairs))
print("INFO: Analysis done in {} pairs".format(len(done)))
print("INFO: {}  had 2 or more LOH reported".format(len(positive)))
