#!/usr/bin/python
# -*- coding: utf-8 -*-

# MAIN program
# Transforms the GISTIC file from TCGA to the same format as ASCAT2, ascatNGS, DNAcopy...
# IMPORTANT: Before running, download OV.focal_score_by_genes.txt from TCGA (https://portal.gdc.cancer.gov/files/14a0766c-6ca4-47bb-ac70-62133c30c1c5)

import os
import shlex
import sqlite3
import subprocess

gistic = "OV.focal_score_by_genes.txt"
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")

if os.path.isfile(gistic) :
    tab = [] # File content
    genes = {} # Equivalence genes: gene_coordinates
    submitters = {} # Equivalence uuid : submitter_id
    order = [] # uuid order according to OV.focal_score_by_genes.txt header
    print("INFO: Reading GISTIC file")
    with open(gistic, "r") as fi :
        for l in fi :
            aux = l.strip().split("\t")
            gene = aux[0].split(".")[0]
            if gene.startswith("ENS") :
                tab.append(aux)
                genes[gene] = {}
            else :
                order = aux[3:]
    print("INFO: Extracting genomic data from UCSC database")
    cmd = "mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -P 3306 -sN -D hg38 -e \"select chrom, chromStart, chromEnd, transcript, protein from knownCanonical\""
    args = shlex.split(cmd)
    proc = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    out, error = proc.communicate()
    auxdb = out.decode().strip().split("\n")
    print("INFO: Adding the coordinates to ENSEMBL id dict")
    # Create a dict with the equivalence ENSG_id -> chromosomal coordinates
    for a in auxdb :
        chr, start, end, transcript, auxgene = a.split("\t")
        g = auxgene.split(".")[0]
        if g in genes.keys() :
            genes[g] = {"chr" : chr, "start" : start, "end" : end}
    # Get the equivalence uuid - submitter_id
    print("INFO: Getting submitters' information from TCGA database")
    with dbcon :
        cur = dbcon.cursor()
        q = cur.execute("SELECT caseId, submitter FROM patient WHERE cancer='OV'")
        cases = q.fetchall()
    for c in cases :
        submitters[c[1]] = {"submitter" : c[0], "content" : "Chromosome\tStart\tEnd\tTotal_CN\tMinor_CN\n"}

    print("INFO: Converting the GISTIC to common format")
    notfound = []
    unknowgene = []
    for t in tab :
        currentgene = t[0].split(".")[0]
        coords = genes[currentgene]
        it = 3
        for o in order :
            aux = t[it]
            it += 1
            tcn = -1
            lcn = -1
            if aux == '-1' :
                tcn = 1
                lcn = 0
            elif aux == '0' :
                tcn = 2
                lcn = 1
            elif aux == '1' :
                tcn = 3
                lcn = 1
            else :
                print("WARNING: {} not a valid value".format(aux))
            if o in submitters.keys() :
                if coords != {} :
                    submitters[o]["content"] += "{chr}\t{sta}\t{end}\t{tcn}\t{lcn}\n".format(chr = coords["chr"], sta = coords["start"], end = coords["end"], tcn = tcn, lcn = lcn)
                elif currentgene not in unknowgene :
                    unknowgene.append(currentgene)
            else :
                if o not in notfound :
                    notfound.append(o)

    print("{} submitters not found".format(len(notfound)))
else :
    print("ERROR: Cannot find GISTIC file in {}".format(gistic))
