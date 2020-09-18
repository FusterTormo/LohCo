#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Check the output by 5 LOH tools (ASCAT2, DNAcopy, FACETS, ascatNGS, and Sequenza) in the genes of interest (BRCA1 and BRCA2)
Assumes that if there is a deleterious mutation in the gene (BRCA1 for example) there will be an LOH in the other gene copy (Deletion or copy number neutral LOH)
Compares how many coincidences each tool reports. A coincidence is when there is the deleterious mutation and the LOH in the gene.
"""

# Libraries
import os
import sqlite3
import sys

import libgetters as lg
import libcomparison as lc
import main1 as lib
import libconstants as cts

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

def getMutation(path, gene) :
    """Classify the worst mutation in the gene

    Given a list of whole exome mutations (path to vcf file), extracts the mutation in the gene, and gets the considered the worst mutation. If this worst mutation is considered deleterious
    returns '+'; if the mutation is considered tolerated returns '-'; in other case the mutation is considered unknown and returns '?'

    Parameters
    ----------
        path : str
            Path to the vcf file that contains all the mutation
        gene : str
            Gene of interest to extract the worst mutation

    Returns
    -------
        str
            * '+' if the worst mutation is considered deleterious ("nonframeshift deletion", "nonframeshift insertion", "frameshift substitution", "frameshift block substitution", "frameshift deletion",
            "frameshift insertion", "stopgain", "splicing")
            * '-' if the worst mutation is considered tolerated ("NA", "synonymous SNV")
            * '?' if the worst mutation is considered unknown ("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift block substitution", "stoploss")
    """
    mut = lib.getWorst(path, gene)
    clas = "-"
    if mut in cts.var_positive :
        clas = "+"
    elif mut in cts.var_negative :
        clas = "-"
    elif mut in cts.var_neutral :
        clas =  "?"
    else :
        clas = "NF"

    return clas

def checkAscat(folder, reg) :
    ascat = "{}/ASCAT2".format(folder)
    cn = ""
    if os.path.isdir(ascat) :
        files = os.listdir(ascat)
        if len(files) == 1 :
            abs = "{}/{}".format(ascat, files[0])
            cn = lib.getLOH(abs, "ascatarray", reg)
        else :
            for f in files :
                abs = "{}/{}".format(ascat, f)
                if cn == "" :
                    cn = lib.getLOH(abs, "ascatarray", reg)
                else :
                    auxCn = lib.getLOH(abs, "ascatarray", reg)
                    if cn != auxCn : # If the ASCAT2 outpus does not output the same aberration, we do not include the result
                        cn =  ""
                        break
    return cn


def doTest(gene, region) :
    output = "ID1\tID2\tgerm\tsom\tASCAT\tFACETS\tascatNGS\tSequenza\n"
    cont = 0
    # Get the OV cases
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
        cases = q.fetchall()

    print("INFO: {} cases found".format(len(cases)))
    # Get bam information from each case and to do the
    for c in cases :
        cont += 1
        if cont % 100 == 0 :
            print("INFO: {} cases checked".format(cont))

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
                workindir = "{wd}/{sub}".format(wd = wd, sub = c[0])
                analysisdir = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]"
                platypust = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
                platypusc = "{}/platypusGerm/platypus.hg38_multianno.txt".format(cf)
                # Get the information regarding the worst variant in the gene selected found in platypus variant calling
                mut_cn = getMutation(platypusc, gene)
                mut_sm = getMutation(platypust, gene)
                # Get the copy number output from ASCAT2
                asc = checkAscat(workindir, region)
                # Get the copy number output from FACETS
                path = "{wd}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = workindir, folder = analysisdir)
                fac = lib.getLOH(path, "facets", region)
                # Get the copy number output from ascatNGS
                if fac == "Not found" :
                    fac = "NF"
                path = lib.findAscatName("{wd}/{folder}_ASCAT/".format(wd = workindir, folder = analysisdir))
                ngs = lib.getLOH(path, "ascatngs", region)
                if ngs == "Not found" :
                    ngs = "NF"
                # Get the copy number output from Sequenza
                path = "{wd}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysisdir, case = c[0], wd = workindir)
                seq = lib.getLOH(path, "sequenza", region)
                if fac == "Not found" :
                    seq = "NF"
                output += "{id1}\t{id2}\t{mtcn}\t{mtsm}\t{cnasc}\t{cnfac}\t{cnngs}\t{cnseq}\n".format(id1 = tm[0], id2 = cn[0], mtcn = mut_cn, mtsm = mut_sm, cnasc = asc, cnfac = fac, cnngs= ngs, cnseq = seq)
    return output

def main() :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    print("INFO: Checking BRCA1")
    txt = doTest("BRCA1", brca1)
    with open("brca1_ascat_facets_ascatngs_sequenza.tsv", "w") as fi :
        fi.write(txt)
    print("INFO: Checking BRCA2")
    txt = doTest("BRCA2", brca2)
    with open("brca2_ascat_facets_ascatngs_sequenza.tsv", "w") as fi :
        fi.write(txt)


if __name__ == "__main__" :
    main()
