#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Check the output by 5 LOH tools (ASCAT2, DNAcopy, FACETS, ascatNGS, and Sequenza) in the genes of interest (BRCA1 and BRCA2)
Assumes that if there is a deleterious mutation in the gene (BRCA1 for example) there will be an LOH in the other gene copy (Deletion or copy number neutral LOH)
Compares how many coincidences each tool reports. A coincidence is when there is the deleterious mutation and the LOH in the gene.
"""

# Libraries
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
        print("WARNING: Unknown mutation found {}".format(mut))

    return clas

def checkAscat(folder, reg) :
    ascat = "{}/ASCAT2".format(folder)
    files = os.listdir(ascat)
    cn = ""
    if len(files) == 1 :
        abs = "{}/{}".format(ascat, f)
        cn = lib.getLOH(abs, "ascatarray", reg)
    else :
        for f in files :
            abs = "{}/{}".format(ascat, f)
            if cn == "" :
                cn = lib.getLOH(abs, "ascatarray", reg)
            else :
                auxCn = lib.getLOH(abs, "ascatarray", reg)
                if cn != auxCn : # If the ASCAT2 outpus does not output the same aberration, we do not include the result
                    cn =  "?"
                    break
    return cn


def doTest(gene, region) :
    # Counters
    total_pos = {"total" : 0, "ascat2" : 0, "facets" : 0, "ascatNGS" : 0, "sequenza" : 0}
    positives = {"ascat2" : 0, "facets" : 0, "ascatNGS" : 0, "sequenza" : 0}
    # Get the OV cases
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
        cases = q.fetchall()

    # Get bam information from each case and to do the
    for c in cases :
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
                platypust = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
                platypusc = "{}/platypusGerm/platypus.hg38_multianno.txt".format(cf)
                # Get the information regarding the worst variant in the gene selected found in platypus variant calling
                mut = getMutation(platypusc, gene)
                if mut == "+" :
                    total_pos["total"] += 1
                    aux = checkAscat(workindir, region)
                    print("{} - {}".format(mut, aux))
                    break
def main() :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    doTest("BRCA1", brca1)

if __name__ == "__main__" :
    main()
