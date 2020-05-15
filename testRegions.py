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

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

def doTest(gene, region) :
    # Gene coordinates, extracted from biogps
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
                # Get the information regarding the worst variant in the gene selected found in platypus variant calling
                vpt1 = lib.getWorst(platypust, gene)
                vpc1 = lib.getWorst(platypusc, gene)
                # Get the LOH information from the different programs
                analysis = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]"
                facets = "{wd}/{sub}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = wd, sub = c[0], folder = analysis)
                lohF = lib.getLOH(facets, "facets", region)
                ascat = lib.findAscatName("{wd}/{case}/{folder}_ASCAT/".format(wd = wd, case = c[0], folder = analysis))
                lohA = lib.getLOH(ascat, "ascatngs", region)
                sequenza = "{wd}/{case}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysis, case = c[0], wd = wd)
                lohS = lib.getLOH(sequenza, "sequenza", region)
                # Check the group for the gene selected
                if vpc1 in positive :
                    # LOH found in positive cases means +1 to tool score
                    cases_positive.append(analysis)
                    if lohF == "L" :
                        scoreF += 1
                    if lohA == "L" :
                        scoreA += 1
                    if lohS == "L" :
                        scoreS += 1
                elif vpc1 in negative :
                    # LOH found in negative cases means -1 to tool score
                    cases_negative.append(analysis)
                    if lohF == "L" :
                        scoreF -= 1
                    if lohA == "L" :
                        scoreA -= 1
                    if lohS == "L" :
                        scoreS -= 1
                else :
                    # If the variant found is nonsynonymous SNV we cannot classify the case. So no score is made
                    cases_neutral.append(analysis)
    results = "INFO: Final score for {}\n".format(gene)
    results += "\tPositive cases: {}\n\tNegative cases: {}\n\tNeutral cases: {}\n".format(len(cases_positive), len(cases_negative), len(cases_neutral))
    results += "\n\tFACETS score: {}\n\tascatNGS score: {}\n\tSequenza score: {}\n".format(scoreF, scoreA, scoreS)
    return results

def main() :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    palb2 = ["16", 23603160, 23641310]
    atm = ["11", 108222484, 108369102]
    print("INFO: Checking BRCA1")
    test1 = doTest("BRCA1", brca1)
    print("INFO: Checking BRCA2")
    test2 = doTest("BRCA2", brca2)
    print("INFO: Checking ATM")
    test3 = doTest("ATM", atm)
    print("INFO: Checking PALB2")
    test4 = doTest("PALB2", palb2)
    print(test1)
    print(test2)
    print(test3)
    print(test4)

main()
