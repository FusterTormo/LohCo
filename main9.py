#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import subprocess
import sys

# Constants
brca1 = ["17", 43044295, 43170245] # Gene of interest to search LOH
cFolder = "fsupek_cancer2"
cancer = "OV" # Cancer repository to search
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")

wd = "/g/strcombio/{cancer_path}/TCGA_bam/{c}".format(c = cancer, cancer_path = cFolder)

# Stats
pairs = 0
done = 0


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
        for cn in controls :
            pairs += 1
            loh = [] # Store the LOH reported by the tools
            prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
            # Get FACETS output file
            folder = "{wd}/{sub}/{pre}_FACETS".format(wd = wd, sub = c[0], pre = prefix)
            file = "{fld}/facets_comp_cncf.tsv".format(fld = folderº)
            if os.path.isfile(folder) :
                # Check LOH in the gene, add the output to a list
                loh.append("L")
            # Get Sequenza output file
            folder = "{wd}/{sub}/{pre}_Sequenza".format(wd = wd, sub = c[0], pre = prefix)
            file = "{fld}/{case}_segments.txt".format(fld = folder, case = c[0])
            if os.path.isfile(file) :
                # Check LOH in the gene, add th output to a list
                loh.append("L")
            # Get ASCAT2 output file
            folder = "{wd}/{sub}/ASCAT2".format(wd = wd, sub = c[0])
            if os.path.isfolder(folder) :
                # Check LOH in the gene, add the output to a list
                loh.append("L")
            if len(loh) == 3 :
                done += 1


        # Count the number of LOH found in the patient
        # https://appdividend.com/2021/02/08/python-list-count-how-to-count-elements-in-python-list/
# Check if FACETS/Sequenza/ASCAT2 have reported LOH
print("INFO: Pairs tumor-control: {}".format(pairs))
print("INFO: Analysis done in {} pairs".format(done))
