#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
TEST TOOL MATCHES
    For each sample, calculates the percent of similarity. Similarity is defined as the number of regions where two (or three) tools report the same
    copy number divided by the total of regions that the case has. Value should be similar to accuracy in the confusion matrix
"""

import os
import sqlite3

import main1 as mm
import libcomparison as lc

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

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
            # Get the absolute path the and the prefix for the tool output
            tf = "{wd}/{sub}/{tm}_VS_{cn}".format(wd = wd, sub = c[0], tm = tm[0].split("-")[0], cn = cn[0].split("-")[0])
            facets = "{}_FACETS/facets_comp_cncf.tsv".format(tf)
            ascat = mm.findAscatName("{}_ASCAT/".format(tf))
            sequenza = "{}_Sequenza/{}_segments.txt".format(tf, c[0])
            if os.path.isfile(facets) :
                outf = lc.convert2region(facets)
            if os.path.isfile(ascat) :
                outa = lc.convert2region(ascat)
            if os.path.isfile(sequenza) :
                outs = lc.convert2region(sequenza)

            if os.path.isfile(facets) and os.path.isfile(ascat) :
                regs = lc.getFragments(outf, outa)
                dc = compi.doComparison(regs, outf, outa)
                print(dc)
            if os.path.isfile(facets) and os.path.isfile(sequenza) :
                regs = lc.getFragments(outf, outs)
                dc = compi.doComparison(regs, outf, outs)
                print(dc)
            if os.path.isfile(ascat) and os.path.isfile(sequenza) :
                regs = lc.getFragments(outa, outs)
                dc = compi.doComparison(regs, outa, outs)
                print(dc)

            
            sys.exit()
