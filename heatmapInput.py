#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Prepares files to be readed by pcaMain6, which will plot heatmap and PCA.
Reads output files from main6
Extracts the best Jaccard Index for each submitter. To extract other column, change variable in line 22 (col)
Stores the data (using ASCAT2 and DNAcopy as truth set) in text files. With the columns:
    Submitter   ASCAT2  DNAcopy FACETS  ascatNGS    Sequenza    PURPLE
"""

import os
import sqlite3
import subprocess

def getJCC(sub, fil) :
    cmd = "grep {} {}".format(sub, fil)
    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = pr.communicate()
    out = std.decode().strip().split("\n")
    header = True
    col = 10 # Column number for Jaccard index
    dat = -1
    for l in out :
        aux = l.split("\t")
        try :
            if l != "" :
                dat = round(float(aux[col]), 4)
        except ValueError :
            pass

    if dat == -1 :
        dat = "NA"

    return dat

"""Main program"""

dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
ascatline = "submitter\tascat2\tdnacopy\tfacets\tascatngs\tsequenza\tpurple\n"
arrayline = "submitter\tascat2\tdnacopy\tfacets\tascatngs\tsequenza\tpurple\n"
# Get submitters list
with dbcon :
      cur = dbcon.cursor()
      q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
      cases = q.fetchall()

for c in cases :
    # For each submitter, extract the best JCC
    ascat = getJCC(c[0], "main6/arrayVSascat2.tsv")
    dnacopy = getJCC(c[0], "main6/arrayVSarray.tsv")
    facets = getJCC(c[0], "main6/arrayVSfacets.tsv")
    ascatngs = getJCC(c[0], "main6/arrayVSascatNGS.tsv")
    sequenza = getJCC(c[0], "main6/arrayVSsequenza.tsv")
    purple = getJCC(c[0], "main6/arrayVSpurple.tsv")
    arrayline += "{sub}\t{a}\t{d}\t{f}\t{n}\t{s}\t{p}\n".format(sub = c[0], a = ascat, d = dnacopy, f = facets, n = ascatngs, s = sequenza, p = purple)

    ascat = getJCC(c[0], "main6/ascat2VSascat2.tsv")
    dnacopy = getJCC(c[0], "main6/ascat2VSarray.tsv")
    facets = getJCC(c[0], "main6/ascat2VSfacets.tsv")
    ascatngs = getJCC(c[0], "main6/ascat2VSascatNGS.tsv")
    sequenza = getJCC(c[0], "main6/ascat2VSsequenza.tsv")
    purple = getJCC(c[0], "main6/ascat2VSpurple.tsv")
    ascatline += "{sub}\t{a}\t{d}\t{f}\t{n}\t{s}\t{p}\n".format(sub = c[0], a = ascat, d = dnacopy, f = facets, n = ascatngs, s = sequenza, p = purple)

with open("main6/dnacopy_heatmap.tsv", "w") as fi :
    fi.write(arrayline)

print("INFO: Data written as main6/dnacopy_heatmap.tsv")

with open("main6/ascat2_heatmap.tsv", "w") as fi :
    fi.write(ascatline)

print("INFO: Data written as main6/ascat2_heatmap.tsv")
print("INFO: Run pcaMain6.R to get PCA and heatmap from this data")
