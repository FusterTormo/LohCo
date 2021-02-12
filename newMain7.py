#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import sys

import main1 as lib

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

# Get submitters list
with dbcon :
      cur = dbcon.cursor()
      q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
      cases = q.fetchall()

for c in cases :
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
        tumors = q.fetchall()
        q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
        controls = q.fetchall()

    # Get the pair tumor/control with more files
    submitter = {"tumor" : "", "control" : "", "vcfFiles" : [], "lohFiles" : []}
    for tm in tumors :
        for cn in controls :
            auxDc = {"tumor" : tm[0], "control" : cn[0], "vcfFiles" : [], "lohFiles" : []}
            tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
            cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])
            workindir = "{wd}/{sub}".format(wd = wd, sub = c[0])
            analysisdir = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]""
            vcf = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
            vcc = "{}/platypusGerm/platypus.hg38_multianno.txt".format(cf)
            # Check both variant calling files exist
            if os.path.isfile(vcf) :
                auxDc["vcfFiles"].append(vcf)
            if os.path.isfile(vcc) :
                auxDc["vcfFiles"].append(vcc)
            folder = "{}/ASCAT2".format(workindir)
            if os.path.isdir(folder) and len(os.listdir(folder)) > 0:
                auxDc["lohFiles"].append(folder)
            facets = "{wd}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = workindir, folder = analysisdir)
            if os.path.isfile(facets) :
                auxDc["lohFiles"].append(facets)
            ascatngs = lib.findAscatName("{wd}/{folder}_ASCAT/".format(wd = workindir, folder = analysisdir))
            if ascatngs != "Not found" :
                auxDc["lohFiles"].append(ascatngs)
            sequenza = "{wd}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysisdir, case = c[0], wd = workindir)
            if os.path.isfile(sequenza) :
                auxDc["lohFiles"].append(sequenza)
            if len(auxDc["vcfFiles"]) > len(submitter["vcfFiles"]) and len(auxDc["lohFiles"]) > len(submitter["vcfFiles"]) :
                submitter = auxDc.copy()

    print(c)
    print(submitter)
