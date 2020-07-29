#!/usr/bin/python
# -*- coding: utf-8 -*-

"""MAIN PROGRAM
TEST TOOL MATCHES
    Calculates the similarity between two tools (FACETS, ascatNGS, and Sequenza)
    When two regions have the same aberration reported by both tools, this is annotated as a coincidence.
    Two different similarities are calculated
        regSimilarity annotates the regions in common that are similar
        baseSimilarity annotates the length of the regions that are annotated as similar

    Additionally, calculates the MCC, and the Jaccard index in all the aberrations
    Finally it calculates the percentage of coincidences by inspecting the regions in common
    The output is printed in separated files (called facetsVSascatngs.tsv, facetsVSsequenza.tsv, and ascatVSsequenza.tsv) one for each comparison
"""

import os
import sqlite3

import libcomparison as lc
import libconstants as cte
import libgetters as lg
import libstatistics as ls
import main1 as mm

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

def main () :
    fvaFi = "facetsVSascatngs.tsv"
    fvsFi = "facetsVSsequenza.tsv"
    avsFi = "ascatVSsequenza.tsv"
    # Write the output files' header
    with open(fvaFi, "w") as fi :
        fi.write("Case\tregSim\tbaseSim\tMCCA\tMCCN\tMCCL\tMCCD\tjcc\tFpurity\tApurity\tFploidy\tAploidy\n")
    with open(fvsFi, "w") as fi :
        fi.write("Case\tregSim\tbaseSim\tMCCA\tMCCN\tMCCL\tMCCD\tjcc\tFpurity\tSpurity\tFploidy\tSploidy\n")
    with open(avsFi, "w") as fi :
        fi.write("Case\tregSim\tbaseSim\tMCCA\tMCCN\tMCCL\tMCCD\tjcc\tApurity\tSpurity\tAploidy\tSploidy\n")

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
                fva.append(tf.split("/")[-1])
                fvs.append(tf.split("/")[-1])
                avs.append(tf.split("/")[-1])
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
                    sts = ls.doContingency(c2) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2)
                    fva.append(ls.regSimilarity(regs, outf, outa))
                    fva.append(ls.baseSimilarity(regs, outf, outa))
                    for ab in cte.aberrations :
                        fva.append(sts[ab]["MCC"])
                    fva.append(jcc)
                    fva.append(outf["purity"])
                    fva.append(outs["purity"])
                    fva.append(outf["ploidy"])
                    fva.append(outs["ploidy"])
                else :
                    fva.append("NA")
                    fva.append("NA")
                    for ab in cte.aberrations :
                        fva.append("NA")
                    for ab in cte.aberrations :
                        fva.append("NA")
                    fva.append("NA")
                    fva.append("NA")
                    fva.append("NA")
                    fva.append("NA")
                # Compare FACETS VS Sequenza
                if os.path.isfile(facets) and os.path.isfile(sequenza) :
                    regs = lc.getFragments(outf, outs)
                    c1 = lc.doComparison(regs, outf, outs)
                    c2 = lc.doComparison2(regs, outf, outs)
                    sts = ls.doContingency(c2) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2) # Get the Jaccard index for all the aberrations
                    fvs.append(ls.regSimilarity(regs, outf, outs))
                    fvs.append(ls.baseSimilarity(regs, outf, outs))
                    for ab in cte.aberrations :
                        fvs.append(sts[ab]["MCC"])
                    fvs.append(jcc)
                    fvs.append(outf["purity"])
                    fvs.append(outs["purity"])
                    fvs.append(outf["ploidy"])
                    fvs.append(outs["ploidy"])
                else :
                    fvs.append("NA")
                    fvs.append("NA")
                    for ab in cte.aberrations :
                        fvs.append("NA")
                    for ab in cte.aberrations :
                        fvs.append("NA")
                    fvs.append("NA")
                    fvs.append("NA")
                    fvs.append("NA")
                    fvs.append("NA")
                # Compare ascatNGS VS Sequenza
                if os.path.isfile(ascat) and os.path.isfile(sequenza) :
                    regs = lc.getFragments(outa, outs)
                    c1 = lc.doComparison(regs, outa, outs)
                    c2 = lc.doComparison2(regs, outa, outs)
                    sts = ls.doContingency(c2) # Get the MCC for all the aberrations
                    jcc = ls.jaccardIndex(c2) # Get the Jaccard index for all the aberrations
                    avs.append(ls.regSimilarity(regs, outa, outs))
                    avs.append(ls.baseSimilarity(regs, outa, outs))
                    for ab in cte.aberrations :
                        avs.append(sts[ab]["MCC"])
                    avs.append(jcc)
                    avs.append(outf["purity"])
                    avs.append(outs["purity"])
                    avs.append(outf["ploidy"])
                    avs.append(outs["ploidy"])
                else :
                    avs.append("NA")
                    avs.append("NA")
                    for ab in cte.aberrations :
                        avs.append("NA")
                    for ab in cte.aberrations :
                        avs.append("NA")
                    avs.append("NA")
                    avs.append("NA")
                    avs.append("NA")
                    avs.append("NA")
                # Write the output in the corresponding files for each comparison
                with open(fvaFi, "a") as fi :
                    fi.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(fva[0], fva[1], fva[2], fva[3], fva[4], fva[5], fva[6], fva[7], fva[8], fva[9], fva[10], fva[11]))
                with open(fvsFi, "a") as fi :
                    fi.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(fvs[0], fvs[1], fvs[2], fvs[3], fvs[4], fvs[5], fvs[6], fvs[7], fvs[8], fvs[9], fvs[10], fvs[11]))
                with open(avsFi, "a") as fi :
                    fi.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(avs[0], avs[1], avs[2], avs[3], avs[4], avs[5], avs[6], avs[7], avs[8], avs[9], avs[10], avs[11]))


if __name__ == "__main__" :
    main()
