#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
TEST REGIONS:
    For each case
        Get the most significative variant in BRCA1, BRCA2, ATM, and PALB2
        Add a classification, according to the variant (+, -, ?)
        Get the CNV/LOH for each region from FACETS, ascatNGS, and Sequenza.
        Score +1 in the positive score if the program has detected LOH in a deleterious mutation
        Score +1 in the negative score if the program has detected LOH in a non-deleterious mutation
"""

import sqlite3
import sys

import libgetters as lg
import libcomparison as lc
import main1 as lib

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"
# Gene coordinates, extracted from biogps
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
palb2 = ["16", 23603160, 23641310]
atm = ["11", 108222484, 108369102]
# Variants classifier
positive = ["nonframeshift substitution", "nonframeshift block substitution", "nonframeshift deletion", "nonframeshift insertion",
    "frameshift substitution", "frameshift block substitution", "frameshift deletion", "frameshift insertion",
    "stopgain"]
negative = ["NA", "synonymous SNV"]

with dbcon :
    cur = dbcon.cursor()
    q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
    cases = q.fetchall()

for c in cases :
    # Recollir la informacio dels bams i el sexe que te el cas registrats
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT uuid FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
        tumors = q.fetchall()
        q = cur.execute("SELECT uuid FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
        controls = q.fetchall()
    for tm in tumors :
        for cn in controls :
            uuid = "{tm}_VS_{cn}".format(tm = tm[0].split("-")[0], cn = cn[0].split("-")[0])
            linea = "{}\t".format(uuid)
            # Get the variant annotation file name
            tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
            cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])
            platypust = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
            platypusc = "{}/platypusGerm/platypus.hg38_multianno.txt".format(cf)
            # Get the FACETS, ascatNGS and sequenza output in REGION format
            ficFa = "{wd}/{sub}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = wd, sub = c[0], folder = analysis)
            if os.path.isfile(ficFa) :
                regFa = lc.convert2region(ficFa, "facets")
            else :
                regFa = "X"
            ficAs = lib.findAscatName("{wd}/{case}/{folder}_ASCAT/".format(wd = wd, case = c[0], folder = analysis))
            if os.path.isfile(ficAs) :
                regAs = lc.convert2region(ficAs, "ascatngs")
            else :
                regAs = "X"
            ficSe = "{wd}/{case}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysis, case = c[0], wd = wd)
            if os.path.isfile(ficSe) :
                regSe = lc.convert2region(ficSe, "sequenza")
            else :
                regSe = "X"
            # Get the information regarding the worst variant in the gene selected found in platypus variant calling
            gene1 = lib.getWorst(platypusc, "BRCA1")
            linea += "{}\t".format(gene1)
            if gene1 in positive :
                linea += "+\t"
            elif gene1 in negative :
                linea += "-\t"
            else :
                linea += "?\t"

            # Get the output from the LOH tools
            if regFa == "X" :
                linea += "-\t"
            else :
                linea += "{}\t".format(lg.getCopyNumber(brca1[1:3], brca1[0], regFa))
            if regAs == "X" :
                linea += "-\t"
            else :
                linea += "{}\t".format(lg.getCopyNumber(brca1[1:3], brca1[0], regAs))
            if regSe == "X" :
                linea += "-\t"
            else :
                linea += "{}\t".format(lg.getCopyNumber(brca1[1:3], brca1[0], regSe))
            print(linea)
            sys.exit()

            # gene2 = lib.getWorst(platypusc, "BRCA2")
            # gene3 = lib.getWorst(platypusc, "ATM")
            # gene4 = lib.getWorst(platypusc, "PALB2")
            #
            # ascat = lib.findAscatName("{wd}/{case}/{folder}_ASCAT/".format(wd = wd, case = c[0], folder = analysis))
            # lohA = lib.getLOH(ascat, "ascatngs", region)
            # sequenza = "{wd}/{case}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysis, case = c[0], wd = wd)
            # lohS = lib.getLOH(sequenza, "sequenza", region)
