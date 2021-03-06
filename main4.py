#!/usr/bin/python
# -*- coding: utf-8 -*-

"""MAIN PROGRAM
TEST REGIONS:
    For each gene
        Get the most significative variant in this gene
        Test if the sample can be included in the positive, negative or unknown group.
        Get the CNV/LOH for this region from FACETS, ascatNGS, and Sequenza.
        Count LOH in each group
            score[] annotates the LOH in positive cases
            score[]2 annotates the LOH in negative cases
            score[]3 annotates the LOH in unknown cases
        [] is the first letter of the tool (i.e. scoreA is ascatNGS score for positive cases, scoreS3 is Sequenza score for unknown cases)
    Data is created for all BRCA1, BRCA2, ATM, and PALB2
    Scores will be printed, and stored in score.txt file
"""

import sqlite3
import sys

import libgetters as lg
import libcomparison as lc
import main1 as lib
import libconstants as cts

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

def doTest(gene, region) :
    """Find the output from FACETS, ascatNGS, and Sequenza to get if there is a LOH in the gene passed as parameter

    Classifies the case as:
        Positive if there is a deleterious variant in the gene (nonframeshift, frameshift, stopgain)
        Negative if there is not a deleterious variant in the gen
        Neutral if the worst variant is nonsynonymous SNV, splicing, or stoploss

    Parameters
    ----------
        gene : str
            Gene name where to find the worst variant and the copy number detected by FACETS, ascatNGS and Sequenza
        region : list
            Coordinates of the gene. The format of the list must be [chrNumber, start, end]

    Returns
    -------
        str
            Summary statistics in a fancy manner
    """
    # Gene coordinates, extracted from biogps
    positive = cts.var_positive
    negative = cts.var_negative
    # Variants considered unknown, so they will be considered as neutral: nonsynonymous SNV, splicing, stoploss
    positiveF = []
    negativeF = []
    neutralF = []
    positiveA = []
    negativeA = []
    neutralA = []
    positiveS = []
    negativeS = []
    neutralS = []
    scoreF = 0
    scoreF2 = 0
    scoreF3 = 0
    scoreA = 0
    scoreA2 = 0
    scoreA3 = 0
    scoreS = 0
    scoreS2 = 0
    scoreS3 = 0

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
                    # Add the case to the corresponding totals, depending if the tool has reported anything
                    if lohF != "Not found" :
                        positiveF.append(analysis)
                    if lohA != "Not found" :
                        positiveA.append(analysis)
                    if lohS != "Not found" :
                        positiveS.append(analysis)

                    # Deleterious variant found. LOH should be present
                    if lohF == "L" :
                        scoreF += 1
                    if lohA == "L" :
                        scoreA += 1
                    if lohS == "L" :
                        scoreS += 1
                elif vpc1 in negative :
                    if lohF != "Not found" :
                        negativeF.append(analysis)
                    if lohA != "Not found" :
                        negativeA.append(analysis)
                    if lohS != "Not found" :
                        negativeS.append(analysis)

                    # Tolerated variant found. LOH should not be present
                    if lohF == "L" :
                        scoreF2 += 1
                    if lohA == "L" :
                        scoreA2 += 1
                    if lohS == "L" :
                        scoreS2 += 1
                else :
                    # If the variant found is nonsynonymous SNV we cannot classify the case. So no score is made
                    if lohF != "Not found" :
                        neutralF.append(analysis)
                    if lohA != "Not found" :
                        neutralA.append(analysis)
                    if lohS != "Not found" :
                        neutralS.append(analysis)

                    # Unknown variant
                    if lohF == "L" :
                        scoreF3 += 1
                    if lohA == "L" :
                        scoreA3 += 1
                    if lohS == "L" :
                        scoreS3 += 1

    results = "INFO: Final score for {}\n".format(gene)
    results += "\tFACETS\n\t\tPositive: {} - Detected: {}\n\t\tNegative: {} - Detected: {}\n\t\tUnknown: {} - Detected {}\n".format(len(positiveF), scoreF, len(negativeF), scoreF2, len(neutralF), scoreF3)
    results += "\tascatNGS\n\t\tPositive: {} - Detected: {}\n\t\tNegative: {} - Detected: {}\n\t\tUnknown: {} - Detected {}\n".format(len(positiveA), scoreA, len(negativeA), scoreA2, len(neutralA), scoreA3)
    results += "\tSequenza\n\t\tPositive: {} - Detected: {}\n\t\tNegative: {} - Detected: {}\n\t\tUnknown: {} - Detected {}\n\n".format(len(positiveS), scoreS, len(negativeS), scoreS2, len(neutralS), scoreS3)
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
    with open("score.txt", "w") as fi :
        fi.write(test1)
        fi.write(test2)
        fi.write(test3)
        fi.write(test4)
    print(test1)
    print(test2)
    print(test3)
    print(test4)
    print("INFO: Information stored in score.txt")

main()
