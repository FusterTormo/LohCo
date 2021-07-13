#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Get, the number of variants reported by Strelka2, and Platypus in a gene defined by the user
    * number of variants found in the gene
    * variant classification
"""

import libconstants as cte

import subprocess
import sqlite3
import os

cancers = {"OV" : "/g/strcombio/fsupek_cancer2/TCGA_bam", "BRCA" : "/g/strcombio/fsupek_cancer3/TCGA_bam"}

def main() :
    dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
    repo = input("Which cancer check the tool?\nOptions: ({}) ".format(" - ".join(cancers.keys())))
    gene = input("Get the variants from which gene? ")
    platy = {"posT" : 0, "posC" : 0, "negT" : 0, "negC" : 0} # Variant classification in Platypus reported variants
    elka = {"posT" : 0, "posC" : 0, "negT" : 0, "negC" : 0}
    pall = {"Tumor" : 0, "Control" : 0} # Total variants detected by Platypus
    sall = {"Tumor" : 0, "Control" : 0}  # Total variants detected by Strelka2
    tumors = 0 # Bam files corresponding to tumor samples
    controls = 0 # Bam files corresponding to control samples
    sample = ""

    # Get the submitter and uuid from the cancer asked
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT s.submitter, uuid, tumor FROM patient p join sample s on s.submitter=p.submitter WHERE cancer='{cancer}'".format(cancer = repo))
        cases = q.fetchall()

    print("INFO: {} bam files found".format(len(cases)))
    for c in cases :
        if c[2].find("Tumor") > 0 :
            tumors += 1
            sample = "Tumor"
        elif c[2].find("Normal") > 0 :
            controls += 1
            sample = "Control"
        else :
            sample = "Undefined"

        if cases.index(c) % 100 == 0 :
            print("INFO: {} cases executed: {} tumors, {} controls".format(cases.index(c), tumors, controls))

        ficP = "{wd}/{cancer}/{sub}/{uuid}/platypusGerm/platypus.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        ficS = "{wd}/{cancer}/{sub}/{uuid}/strelkaGerm/results/variants/strelka.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        if os.path.isfile(ficP) : # Check if platypus file exists
            cmd = "grep {gene} {file}".format(gene = gene, file = ficP)
            pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std, err = pr.communicate()
            out = std.decode().strip().split("\n")
            for o in out :
                pall[sample] += 1
                aux = o.split("\t")
                if len(aux) > 1 :
                    if aux[5].find("splicing") > 0 :
                        if sample == "Tumor" :
                            platy["posT"] += 1
                        elif sample == "Control" :
                            platy["posC"] += 1
                    else :
                        if aux[8] in cte.var_positive :
                            if sample == "Tumor" :
                                platy["posT"] += 1
                            elif sample == "Control" :
                                platy["posC"] += 1
                        elif aux[8] in cte.var_negative :
                            if sample == "Tumor" :
                                platy["negT"] += 1
                            elif sample == "Control" :
                                platy["negC"] += 1
        else :
            print("WARNING: {} not found".format(ficP))

        if os.path.isfile(ficS) :
            cmd = "grep {gene} {file}".format(gene = gene, file = ficS)
            pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std, err = pr.communicate()
            out = std.decode().strip().split("\n")
            for o in out :
                sall[sample] += 1
                aux = o.split("\t")
                if len(aux) > 1 :
                    if aux[5].find("splicing") > 0 :
                        if sample == "Tumor" :
                            elka["posT"] += 1
                        elif sample == "Control" :
                            elka["posC"] += 1
                    else :
                        if aux[8] in cte.var_positive :
                            if sample == "Tumor" :
                                elka["posT"] += 1
                            elif sample == "Control" :
                                elka["posC"] += 1
                        elif aux[8] in cte.var_negative :
                            if sample == "Tumor" :
                                elka["negT"] += 1
                            elif sample == "Control" :
                                elka["negC"] += 1
        else :
            print("WARNING: {} not found".format(ficS))

    print("SUMMARY\n-------")
    print("Platypus reported {} variants in {} gene".format(pall, gene))
    print("In TUMOR samples: {} Positive, {} Negative".format(platy["posT"], platy["negT"]))
    print("In CONTROL samples: {} Positve, {} Negative\n".format(platy["posC"], platy["negC"]))
    print("Strelka reported {} variants in {} gene".format(sall, gene))
    print("In TUMOR samples: {} Positive, {} Negative".format(elka["posT"], elka["negT"]))
    print("In CONTROL samples: {} Positve, {} Negative".format(elka["posC"], elka["negC"]))
if __name__ == "__main__" :
    main()
