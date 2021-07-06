#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Get, the number of variants reported by Strelka2, and Platypus in a gene defined by the user
    * number of variants found in the gene
    * variant classification
"""

import sqlite3
import os

cancers = {"OV" : "/g/strcombio/fsupek_cancer2/TCGA_bam", "BRCA" : "/g/strcombio/fsupek_cancer3/TCGA_bam"}

def main() :
    dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
    repo = input("Which cancer check the tool?\nOptions: ({}) ".format(" - ".join(cancers.keys())))
    gene = input("Get the variants from which gene? ")

    # Get the submitter and uuid from the cancer asked
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT s.submitter, uuid FROM patient p join sample s on s.submitter=p.submitter WHERE cancer='{cancer}'".format(cancer = repo))
        cases = q.fetchall()

    for c in cases :
        ficP = "{wd}/{cancer}/{sub}/{uuid}/platypusGerm/platypus.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        ficS = "{wd}/{cancer}/{sub}/{uuid}/strelkaGerm/results/variants/strelka.hg38_multianno.txt".format(wd = cancers[repo], cancer = repo, sub = c[0], uuid = c[1])
        if os.path.isfile(ficP) : # Check if folder exists
            cmd = "grep {gene} {file}".format(gene = gene, file = ficP)
            print(cmd)
        else :
            print("WARNING: {} not found".format(dirP))

if __name__ == "__main__" :
    main()
