#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3

import libcomparison as lc
import libgetters as lg
import libstatistics as ls

def getFACETS(path) :
    """ Return the FACETS output files in the submitter folder passed as parameter """
    facets = []
    folders = os.listdir(path)
    for f in folders :
        if f.endswith("_FACETS") :
            aux = "{abs}/{fa}/facets_comp_cncf.tsv".format(abs = path, fa = f)
            if os.path.isfile(aux) :
                facets.append(aux)
    return facets

def getAscatNGS(path) :
    """Return the ascatNGS output files in the submitter folder passed as parameter """
    ascat = []
    folders = os.listdir(path)
    for f in folders :
        if f.endswith("_ASCAT") :
            aux = "{abs}/{ngs}".format(abs = path, ngs = f)
            folders2 = os.listdir(aux)
            for g in folders2 :
                if g.endswith("copynumber.caveman.csv") :
                    aux2 = "{}/{}".format(aux, g)
                    ascat.append(aux2)
    return ascat

def getSequenza(path) :
    """ Return the Sequenza output files in the submitter folder passed as parameter """
    sequenza = []
    folders = os.listdir(path)
    for f in folders :
        if f.endswith("_Sequenza") :
            aux = "{abs}/{sq}".format(abs = path, sq = f)
            folders2 = os.listdir(aux)
            for g in folders2 :
                if g.endswith("_segments.txt") :
                    aux2 = "{}/{}".format(aux, g)
                    sequenza.append(aux2)
    return sequenza

def compareTools(reg1, reg2) :
    regs = lc.getFragments(reg1, reg2)
    comp = lc.doComparison2(regs, reg1, reg2)
    mat = ls.doContingency(comp)
    jcc = ls.jaccardIndex(comp)
    num = ls.regionNumber(regs)
    base = ls.baseSimilarity(regs, reg1, reg2)
    regions = ls.regSimilarity(regs, reg1, reg2)
    st = "{nr}\t{bs}\t{rs}\t{mcca}\t{mccn}\t{mccl}\t{mccd}\t{jcc}".format(
        nr = num, bs = base, rs = regions, mcca = mat["A"]["MCC"], mccn = mat["N"]["MCC"], mccl = mat["L"]["MCC"], mccd = mat["D"]["MCC"], jcc = jcc)
    return st

def createFile(name) :
    """Creates a new file with the name passed as parameter. Prints the header in the file"""
    with open(name, "w") as fi :
        fi.write("ID1\tID2\tREGIONS\tBASEsim\tREGIONsim\tMCCamp\tMCCnorm\tMCCloh\tMCCdel\tJCC\n")

def main() :
    # Constants
    dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
    cancer = "OV"
    cancerpath = "/g/strcombio/fsupek_cancer2/TCGA_bam/"
    # Get the OV submitters from the database
    with dbcon :
        query = "SELECT submitter FROM patient WHERE cancer='{}' LIMIT 2".format(cancer)
        c = dbcon.cursor()
        x = c.execute(query)
        submitters = x.fetchall()

    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        ascatFolder = "{}/ASCAT2/".format(workindir)
        # Open ASCAT2 folder and get the files available
        if os.path.isdir (ascatFolder) :
            ascatFiles = os.listdir(ascatFolder)
            # Open SNP-array folder and get the files that are in
            arrayFolder = "{}/Array/".format(workindir)
            # Compare SNP-Arrays CNV outputs with ASCAT2
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                print("INFO: Comparing ASCAT2 and Array outputs in {}".format(s))
                for a in ascatFiles :
                    ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                    for b in arrayFiles :
                        arr = lc.convert2region("{}/{}".format(arrayFolder, b), "array")
                        if not os.path.isfile("ascatVSarray.tsv") :
                            createFile("ascatVSarray.tsv")
                        with open("ascatVSarray.tsv", "a") as fi :
                            fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(ascat, arr)))


            # Compare FACETS LOH/CNV outputs with ASCAT2
            print("INFO: Comparing ASCAT2 and FACETS outputs in {}".format(s))
            facetsFiles = getFACETS(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in facetsFiles :
                    f = lc.convert2region(b, "facets", "error")
                    print(compareTools(ascat, f))

            # Compare ascatNGS LOH/CNV outputs with ASCAT2
            print("INFO: Comparing ASCAT2 and ascatNGS outputs in {}".format(s))
            ascatngsFiles = getAscatNGS(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in ascatngsFiles :
                    ngs = lc.convert2region(b, "ascatngs", "error")
                    print(compareTools(ascat, ngs))

            # Compare Sequenza LOH/CNV outputs with ASCAT2
            print("INFO: Comparing ASCAT2 and Sequenza outputs in {}".format(s))
            sequenzaFiles = getSequenza(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in sequenzaFiles :
                    s = lc.convert2region(b, "sequenza", "error")
                    print(compareTools(ascat, s))

if __name__ == "__main__" :
    main()
