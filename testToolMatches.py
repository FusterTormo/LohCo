#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
TEST TOOL MATCHES
    For each sample, calculates the percent of similarity. Similarity is defined as the number of regions where two (or three) tools report the same
    copy number divided by the total of regions that the case has. Value should be similar to accuracy in the confusion matrix
    Additionally, calculates the MCC, and the Jaccard index for all the aberrations
    Finally it calculates the percentage of coincidences by inspecting the regions in common
    The output is printed in separated files one for each comparison
"""

import os
import sqlite3

import main1 as mm
import libcomparison as lc
import libstatistics as ls
import libconstants as cte
import libgetters as lg

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

def percentSimilarity(regs, tool1, tool2) :
    coin = 0
    all = 0
    for k in regs.keys() :
        for r in regs[k] :
            all += 1
            if lg.getCopyNumber(r, k, tool1) == lg.getCopyNumber(r, k, tool2) :
                coin += 1

    percent = 100*float(coin)/float(all)
    return percent


def calculateSimilarity(dc) :
    dividendo = 0
    divisor = 0
    for a in cte.aberrations :
        dividendo += dc[a][a]
        for b in cte.aberrations :
            divisor += dc[a][b]
    similarity = 100*float(dividendo)/float(divisor)
    return similarity

def main () :
    fvaFi = "facetsVSascatngs.tsv"
    fvsFi = "facetsVSsequenza.tsv"
    avsFi = "ascatVSsequenza.tsv"
    # Write the output files' header
    with open(fvaFi, "w") as fi :
        fi.write("Case\tpercent\tACC\tMCCA\tMCCN\tMCCL\tMCCD\tjcca\tjccn\tjccl\tjccd\n")
    with open(fvsFi, "w") as fi :
        fi.write("Case\tpercent\tACC\tMCCA\tMCCN\tMCCL\tMCCD\tjcca\tjccn\tjccl\tjccd\n")
    with open(avsFi, "w") as fi :
        fi.write("Case\tpercent\tACC\tMCCA\tMCCN\tMCCL\tMCCD\tjcca\tjccn\tjccl\tjccd\n")

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
                fva = []
                fvs = []
                avs = []
                # Get the absolute path the and the prefix for the tool output
                tf = "{wd}/{sub}/{tm}_VS_{cn}".format(wd = wd, sub = c[0], tm = tm[0].split("-")[0], cn = cn[0].split("-")[0])
                fva.append(tf)
                fvs.append(tf)
                avs.append(tf)
                facets = "{}_FACETS/facets_comp_cncf.tsv".format(tf)
                ascat = mm.findAscatName("{}_ASCAT/".format(tf))
                sequenza = "{}_Sequenza/{}_segments.txt".format(tf, c[0])
                if os.path.isfile(facets) :
                    outf = lc.convert2region(facets, "facets")
                if os.path.isfile(ascat) :
                    outa = lc.convert2region(ascat, "ascatngs")
                if os.path.isfile(sequenza) :
                    outs = lc.convert2region(sequenza, "sequenza")
                # Compare FACETS vs ascatNGS
                if os.path.isfile(facets) and os.path.isfile(ascat) :
                    regs = lc.getFragments(outf, outa)
                    c1 = lc.doComparison(regs, outf, outa)
                    c2 = lc.doComparison2(regs, outf, outa)
                    sts = ls.doContingency(c1) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2)
                    fva.append(str(calculateSimilarity(c2)))
                    fva.append(str(percentSimilarity(regs, outf, outa)))
                    for ab in cte.aberrations :
                        fva.append(str(sts[ab]["MCC"]))
                    for ab in cte.aberrations :
                        fva.append(str(jcc[ab]))
                else :
                    fva.append("NA")
                    fva.append("NA")
                    for ab in cte.aberrations :
                        fva.append("NA")
                    for ab in cte.aberrations :
                        fva.append("NA")
                # Compare FACETS VS Sequenza
                if os.path.isfile(facets) and os.path.isfile(sequenza) :
                    regs = lc.getFragments(outf, outs)
                    c1 = lc.doComparison(regs, outf, outs)
                    c2 = lc.doComparison2(regs, outf, outs)
                    sts = ls.doContingency(c1) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2) # Get the Jaccard index for all the aberrations
                    fvs.append(str(calculateSimilarity(c2)))
                    fvs.append(str(percentSimilarity(regs, outf, outs)))
                    for ab in cte.aberrations :
                        fvs.append(str(sts[ab]["MCC"]))
                    for ab in cte.aberrations :
                        fvs.append(str(jcc[ab]))
                else :
                    fvs.append("NA")
                    fvs.append("NA")
                    for ab in cte.aberrations :
                        fvs.append("NA")
                    for ab in cte.aberrations :
                        fvs.append("NA")
                # Compare ascatNGS VS Sequenza
                if os.path.isfile(ascat) and os.path.isfile(sequenza) :
                    regs = lc.getFragments(outa, outs)
                    c1 = lc.doComparison(regs, outa, outs)
                    c2 = lc.doComparison2(regs, outa, outs)
                    sts = ls.doContingency(c1) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2) # Get the Jaccard index for all the aberrations
                    fvs.append(str(calculateSimilarity(c2)))
                    fvs.append(str(percentSimilarity(regs, outa, outs)))
                    for ab in cte.aberrations :
                        avs.append(str(sts[ab]["MCC"]))
                    for ab in cte.aberrations :
                        avs.append(str(jcc[ab]))
                else :
                    fvs.append("NA")
                    fvs.append("NA")
                    for ab in cte.aberrations :
                        avs.append("NA")
                    for ab in cte.aberrations :
                        avs.append("NA")
                # Write the output in the corresponding files for each comparison
                with open(fvaFi, "a") as fi :
                    fi.write("\t".join(fva))
                    fi.write("\n")
                with open(fvsFi, "a") as fi :
                    fi.write("\t".join(fvs))
                    fi.write("\n")
                with open(avsFi, "a") as fi :
                    fi.write("\t".join(avs))
                    fi.write("\n")



main()
