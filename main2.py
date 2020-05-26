#!/usr/bin/python
# -*- coding: utf-8 -*-

"""MAIN PROGRAM
TEST REGIONS IN COMMON
    Annotates the regions in common for two tools in each case
    When two regions have the same aberration reported by both tools, this region is annotated as a coincidence.
    Otherwise, the region is annotated as a difference.
    From this list of coincidences/differences, get largest and smallest region
"""

import os
import sqlite3
import sys

import libcomparison as lc
import libgetters as lg
import main1 as mm


def compareRegions(regions, tool1, tool2) :
    coincide = []
    differ = []
    cont = 0
    for k in regs.keys() : # Iterate by chromosome
        for r in regs[k] : # Iterate the regions in the chromosome
            cont += 1
            if lg.getCopyNumber(r, k, tool1) == lg.getCopyNumber(r, k, tool2) :
                coincide.append(r)
            else :
                differ.append(r)
    print("INFO: {} regions found".format(cont))
    small = 9999999999999999999999999999999999999999999999999999
    large = -1
    for r in coincide :
        length = int(r[1]) - int(r[0])
        if length > 0 :
            if length > large :
                large = length
            if length < small :
                small = length
    print("{} regions in common. Smallest: {}. Largest: {}".format(len(coincide), small, large))
    small = 9999999999999999999999999999999999999999999999999999
    large = -1
    for r in differ :
        length = int(r[1]) - int(r[0])
        if length > 0 :
            if length > large :
                large = length
            if length < small :
                small = length
    print("{} regions in common. Smallest: {}. Largest: {}".format(len(differ), small, large))

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

table = []
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
            tf = "{wd}/{sub}/{tm}_VS_{cn}".format(wd = wd, sub = c[0], tm = tm[0].split("-")[0], cn = cn[0].split("-")[0])
            facets = "{}_FACETS/facets_comp_cncf.tsv".format(tf)
            ascat = mm.findAscatName("{}_ASCAT/".format(tf))
            sequenza = "{}_Sequenza/{}_segments.txt".format(tf, c[0])
            if os.path.isfile(facets) :
                outf = lc.convert2region(facets, "facets", "quiet")
            if os.path.isfile(ascat) :
                outa = lc.convert2region(ascat, "ascatngs", "quiet")
            if os.path.isfile(sequenza) :
                outs = lc.convert2region(sequenza, "sequenza", "quiet")
            # Compare FACETS vs ascatNGS
            if os.path.isfile(facets) and os.path.isfile(sequenza) :
                regs = lc.getFragments(outf, outs)
                compareRegions(regs, outf, outs)
                sys.exit()
