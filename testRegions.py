#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
TEST REGIONS:
    For each case
        Get the most significative variant
        Test if the sample can be included in the positive, negative or unknown group.
        Get the CNV/LOH for this region from FACETS, ascatNGS, and Sequenza.
        Score +1 if the program has detected LOH in a deleterious mutation
        Score -1 if the program has detected LOH in a non-deleterious mutation
"""

import sqlite3
import sys

import libgetters as lg
import libcomparison as lc
import main1 as lib

# Gene coordinates, extracted from biogps
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
palb2 = ["16", 23603160, 23641310]
atm = ["11", 108222484, 108369102]
positive = ["nonframeshift substitution", "nonframeshift deletion",  "frameshift deletion", "frameshift insertion", "stopgain"]
negative = ["NA", "synonymous SNV"]
cases_positive = []
cases_negative = []
cases_neutral = []
scoreF = 0
scoreA = 0
scoreS = 0

with dbcon :
    cur = dbcon.cursor()
    q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
    cases = q.fetchall()

for c in cases :
    cont += 1
    # Recollir la informacio dels bams i el sexe que te el cas registrats
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
        tumors = q.fetchall()
        q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
        controls = q.fetchall()
    for tm in tumors :
        for cn in controls :
            # Get the absolute path for the platypus hg38_multianno file in tumor and control
            tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
            cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])
            platypust = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
            platypusc = "{}/platypusGerm/platypus.hg38_multianno.txt".format(cf)
            # Get the information regarding the worst variant in BRCA1 found in platypus variant calling
            vpt1 = lib.getWorst(platypust, "BRCA1")
            vpc1 = lib.getWorst(platypusc, "BRCA1")
            # Get the LOH information from the different programs
            analysis = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]"
            facets = "{wd}/{sub}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = wd, sub = c[0], folder = analysis)
            lohF = lib.getLOH(facets, "facets", brca1)
            ascat = lib.findAscatName("{wd}/{case}/{folder}_ASCAT/".format(wd = wd, case = c[0], folder = analysis))
            lohA = getLOH(ascat, "ascatngs", brca1)
            sequenza = "{wd}/{case}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysis, case = c[0], wd = wd)
            lohS = getLOH(sequenza, "sequenza", brca1)
            print("{}\t{}\t{}\t{}".format(vpc1, lohF, lohA, lohS))
            sys.exit()
