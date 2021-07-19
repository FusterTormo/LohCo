#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import sys

import main1 as lib
import libcomparison as lc
import libstatistics as ls

"""MAIN PROGRAM
    Search all the available pairs tumor-control from an optional cancer passed as parameter (default OV).
    For each pair tumor-control, find the tools that have output LOH report.
    Calculate mean copy number and the percentage of (A)mplification, (L)OH, (D)eletion and (N)ormal copy number.
    Output all the information in a tab-delimited file.
        Columns: submitter, case, fac_meanCN, fac_purity, fac_ploidy, fac_aberration, asc_meanCN, asc_aberration, seq_meanCN, seq_purity, seq_ploidy, seq_aberration, pur_meanCN,
        pur_purity, pur_ploidy, pur_aberration, ngs_meanCN, ngs_purity, ngs_ploidy, ngs_aberration
"""

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

# Functions
def convertToCSV(data) :
    return ";".join([str(data["perA"]),str(data["perL"]),str(data["perD"]),str(data["perN"])])

def main(cancer = "OV") :
    """Main program"""

    wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/{}".format(cancer)
    txt = "submitter\tcase\tfac_meanCN\tfac_purity\tfac_ploidy\tfac_aberration\tasc_meanCN\tasc_aberration\tseq_meanCN\tseq_purity\tseq_ploidy\tseq_aberration\tpur_meanCN\tpur_purity\t"
    txt += "pur_ploidy\tpur_aberration\tngs_meanCN\tngs_purity\tngs_ploidy\tngs_aberration\n"
    na = "NA"
    outputFile = "meanCN.tsv"
    count = 0

    # Get submitters list
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='{}'".format(cancer))
          cases = q.fetchall()

    print("INFO: Analysis done in {} cases".format(len(cases)))
    for c in cases :
        count += 1
        if count % 100 == 0 :
            print("INFO: {} cases done".format(count))

        with dbcon :
            cur = dbcon.cursor()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
            tumors = q.fetchall()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
            controls = q.fetchall()

        for tm in tumors :
            for cn in controls :
                tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
                cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])
                workindir = "{wd}/{sub}".format(wd = wd, sub = c[0])
                analysisdir = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]""
                # From each tool get the purity/ploidy
                rAscat = {"purity" : na, "ploidy" : na}
                # And calculate, using libstatistics, the mean copy number, and the percentage of (A)mplifications, (L)OH, (D)eletion or (N)ormal copy number
                sAscat = {"meanCN" : na, "perA" : na, "perL" : na, "perD" : na, "perN" : na}
                rFacets = {"purity" : na, "ploidy" : na}
                sFacets = {"meanCN" : na, "perA" : na, "perL" : na, "perD" : na, "perN" : na}
                rNgs = {"purity" : na, "ploidy" : na}
                sNgs = {"meanCN" : na, "perA" : na, "perL" : na, "perD" : na, "perN" : na}
                rSequenza = {"purity" : na, "ploidy" : na}
                sSequenza = {"meanCN" : na, "perA" : na, "perL" : na, "perD" : na, "perN" : na}
                rPurple = {"purity" : na, "ploidy" : na}
                sPurple = {"meanCN" : na, "perA" : na, "perL" : na, "perD" : na, "perN" : na}

                folder = "{}/ASCAT2".format(workindir)
                # Collect and calculate all the data
                # From ASCAT2
                if os.path.isdir(folder) and len(os.listdir(folder)) > 0:
                    temp = os.listdir(folder)[0] # TODO: Check all ASCAT files
                    ascat = "{wd}/{fi}".format(wd = folder, fi = temp)
                    rAscat = lc.convert2region(ascat, "ascatarray", "error")
                    sAscat = ls.meanCoverage(rAscat)
                # From FACETS
                facets = "{wd}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = workindir, folder = analysisdir)
                if os.path.isfile(facets) :
                    rFacets = lc.convert2region(facets, "facets", "error")
                    sFacets = ls.meanCoverage(rFacets)
                # From ascatNGS
                ascatngs = lib.findAscatName("{wd}/{folder}_ASCAT/".format(wd = workindir, folder = analysisdir))
                if ascatngs != "Not found" :
                    rNgs = lc.convert2region(ascatngs, "ascatngs", "error")
                    sNgs = ls.meanCoverage(rNgs)
                # From Sequenza
                sequenza = "{wd}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysisdir, case = c[0], wd = workindir)
                if os.path.isfile(sequenza) :
                    rSequenza = lc.convert2region(sequenza, "sequenza", "error")
                    sSequenza = ls.meanCoverage(rSequenza)
                # From PURPLE
                purple = "{wd}/{folder}_PURPLE/TUMOR.purple.cnv.somatic.tsv".format(wd = workindir, folder = analysisdir)
                if os.path.isfile(purple) :
                    rPurple = lc.convert2region(purple, "purple", "error")
                    sPurple = ls.meanCoverage(rPurple)

                # Write the output in RAM
                txt += "{sub}\t{an}\t{fmcn}\t{fpu}\t{fpl}\t{fab}\t{acn}\t{aab}\t{scn}\t{spu}\t{spl}\t{sab}\t{pcn}\t{ppu}\t{ppl}\t{pab}\t{ncn}\t{npu}\t{npl}\t{nab}\n".format(
                sub = c[0], an = analysisdir, fmcn = sFacets["meanCN"], fpu = rFacets["purity"], fpl = rFacets["ploidy"], fab = convertToCSV(sFacets),
                acn = sAscat["meanCN"], aab = convertToCSV(sAscat),
                scn = sSequenza["meanCN"], spu = rSequenza["purity"], spl = rSequenza["ploidy"], sab = convertToCSV(sSequenza),
                pcn = sPurple["meanCN"], ppu = rPurple["purity"], ppl = rPurple["ploidy"], pab = convertToCSV(sPurple),
                ncn = sNgs["meanCN"], npu = rNgs["purity"], npl = rNgs["ploidy"], nab = convertToCSV(sNgs))

    with open(outputFile, "w") as fi :
        fi.write(txt)
    print("INFO: Data stored in {} file".format(outputFile))

if __name__ == "__main__" :
    main("OV")
