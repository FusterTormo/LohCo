#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Get, the number of variants reported by Strelka2, and Platypus in a gene defined by the user
    * number of variants found in the gene
    * variant classification
"""

import libconstants as cte

import sqlite3
import os

cancers = {"OV" : "/g/strcombio/fsupek_cancer2/TCGA_bam", "BRCA" : "/g/strcombio/fsupek_cancer3/TCGA_bam"}

def main() :
    dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
    repo = input("Which cancer check the tool?\nOptions: ({}) ".format(" - ".join(cancers.keys())))
    gene = input("Get the variants from which gene? ")
    pall = 0 # Number of variants reported by Platypus
    ppos = 0
    pneg = 0
    sall = 0 # Number of variants reported by Strelka2
    spos = 0
    sneg = 0

    # Get the submitter and uuid from the cancer asked
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT s.submitter, uuid FROM patient p join sample s on s.submitter=p.submitter WHERE cancer='{cancer}'".format(cancer = repo))
        cases = q.fetchall()

    for c in cases :
        ficP = "{wd}/{cancer}/{sub}/{uuid}/platypusGerm/platypus.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        ficS = "{wd}/{cancer}/{sub}/{uuid}/strelkaGerm/results/variants/strelka.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        if os.path.isfile(ficP) : # Check if platypus file exists
            cmd = "grep {gene} {file}".format(gene = gene, file = ficP)
            pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std, err = pr.communicate()
            out = std.decode().strip().split("\n")
            for o in out :
                pall += 1
                aux = o.split("\t")
                if len(aux) > 1 :
                    if aux[5].find("splicing") > 0 :
                        ppos += 1
                    else :
                        if aux[8] in cte.var_positive :
                            ppos += 1
                        elif aux[8] in cte.var_negative :
                            pneg += 1

        else :
            print("WARNING: {} not found".format(ficP))

if __name__ == "__main__" :
    main()
