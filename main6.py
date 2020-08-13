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
        query = "SELECT submitter FROM patient WHERE cancer='{}'".format(cancer)
        c = dbcon.cursor()
        x = c.execute(query)
        submitters = x.fetchall()

    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} ASCAT".format(s))
        ascatFolder = "{}/ASCAT2/".format(workindir)
        if os.path.isdir (ascatFolder) :
            # Open ASCAT2 folder and get the files available
            ascatFiles = os.listdir(ascatFolder)
            # Compare ASCAT2 with itself
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                if not os.path.isfile("ascat2VSascat2.tsv") :
                    createFile("ascat2VSascat2.tsv")
                with open("ascat2VSascat2.tsv", "a") as fi :
                    fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = a, cmp = compareTools(ascat, ascat)))
            # Open SNP-array folder and get the files that are in
            arrayFolder = "{}/Array/".format(workindir)
            # Compare SNP-Arrays CNV outputs with ASCAT2
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                # print("INFO: Comparing ASCAT2 and Array outputs in {}".format(s))
                for a in ascatFiles :
                    ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                    for b in arrayFiles :
                        arr = lc.convert2region("{}/{}".format(arrayFolder, b), "array")
                        if not os.path.isfile("ascat2VSarray.tsv") :
                            createFile("ascat2VSarray.tsv")
                        with open("ascat2VSarray.tsv", "a") as fi :
                            fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(ascat, arr)))


            # Compare FACETS LOH/CNV outputs with ASCAT2
            facetsFiles = getFACETS(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in facetsFiles :
                    f = lc.convert2region(b, "facets", "error")
                    if not os.path.isfile("ascat2VSfacets.tsv") :
                        createFile("ascat2VSfacets.tsv")
                    with open("ascat2VSfacets.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(ascat, f)))

            # Compare ascatNGS LOH/CNV outputs with ASCAT2
            # print("INFO: Comparing ASCAT2 and ascatNGS outputs in {}".format(s))
            ascatngsFiles = getAscatNGS(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in ascatngsFiles :
                    ngs = lc.convert2region(b, "ascatngs", "error")
                    if not os.path.isfile("ascat2VSascatNGS.tsv") :
                        createFile("ascat2VSascatNGS.tsv")
                    with open("ascat2VSascatNGS.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(ascat, ngs)))

            # Compare Sequenza LOH/CNV outputs with ASCAT2
            # print("INFO: Comparing ASCAT2 and Sequenza outputs in {}".format(s))
            sequenzaFiles = getSequenza(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in sequenzaFiles :
                    s = lc.convert2region(b, "sequenza", "error")
                    if not os.path.isfile("ascat2VSsequenza.tsv") :
                        createFile("ascat2VSsequenza.tsv")
                    with open("ascat2VSsequenza.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(ascat, s)))

    # Move the output data to a new folder
    os.mkdir("main6")
    os.rename("ascat2VSascat2.tsv", "main6/ascat2VSascat2.tsv")
    os.rename("ascat2VSarray.tsv", "main6/ascat2VSarray.tsv")
    os.rename("ascat2VSfacets.tsv", "main6/ascat2VSfacets.tsv")
    os.rename("ascat2VSascatNGS.tsv", "main6/ascat2VSascatNGS.tsv")
    os.rename("ascat2VSsequenza.tsv", "main6/ascat2VSsequenza.tsv")

    # Repeat the analysis, but using Arrays as True set
    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} arrays".format(s))
        # Open SNP-Array folder and get the files available
        arrayFolder = "{}/Array/".format(workindir)
        if os.path.isdir (arrayFolder) :
            arrayFiles = os.listdir(arrayFolder)
            # Compare arrays with itself
            for a in arrayFiles :
                arr = lc.conver2region("{}/{}".format(arrayFolder, a), "array")
                if not os.path.isfile("arrayVSarray.tsv") :
                    createFile("arrayVSarray.tsv")
                with open("arrayVSarray.tsv", "a") as fi :
                    fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = a, cmp = compareTools(arr, arr)))
            # Open ASCAT2 folder to get the files that are in
            ascatFolder = "{}/ASCAT2/".format(workindir)
            # Compare ASCAT2 outputs with Arrays
            if os.path.isdir(ascatFolder) :
                # print("INFO: Comparing ASCAT2 and Array outputs in {}".format(s))
                ascatFiles = os.listdir(ascatFolder)
                for a in arrayFiles :
                    arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                    for b in ascatFiles :
                        ascat = lc.convert2region("{}/{}".format(ascatFolder, b), "ascatarray")
                        if not os.path.isfile("arrayVSascat2.tsv") :
                            createFile("arrayVSascat2.tsv")
                        with open("arrayVSascat2.tsv", "a") as fi :
                            fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(arr, ascat)))


            # Compare FACETS LOH/CNV outputs with SNP-Array
            # print("INFO: Comparing ASCAT2 and FACETS outputs in {}".format(s))
            facetsFiles = getFACETS(workindir)
            for a in arrayFiles :
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                for b in facetsFiles :
                    f = lc.convert2region(b, "facets", "error")
                    if not os.path.isfile("arrayVSfacets.tsv") :
                        createFile("arrayVSfacets.tsv")
                    with open("arrayVSfacets.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(arr, f)))

            # Compare ascatNGS LOH/CNV outputs with SNP-Array
            # print("INFO: Comparing ASCAT2 and ascatNGS outputs in {}".format(s))
            ascatngsFiles = getAscatNGS(workindir)
            for a in arrayFiles :
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                for b in ascatngsFiles :
                    ngs = lc.convert2region(b, "ascatngs", "error")
                    if not os.path.isfile("arrayVSascatNGS.tsv") :
                        createFile("arrayVSascatNGS.tsv")
                    with open("arrayVSascatNGS.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(arr, ngs)))

            # Compare Sequenza LOH/CNV outputs with ASCAT2
            # print("INFO: Comparing ASCAT2 and Sequenza outputs in {}".format(s))
            sequenzaFiles = getSequenza(workindir)
            for a in arrayFiles :
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                for b in sequenzaFiles :
                    s = lc.convert2region(b, "sequenza", "error")
                    if not os.path.isfile("arrayVSsequenza.tsv") :
                        createFile("arrayVSsequenza.tsv")
                    with open("arrayVSsequenza.tsv", "a") as fi :
                        fi.write("{id1}\t{id2}\t{cmp}\n".format(id1 = a, id2 = b, cmp = compareTools(arr, s)))

    os.rename("arrayVSarray.tsv", "main6/arrayVSarray.tsv")
    os.rename("arrayVSascat2.tsv", "main6/arrayVSascat2.tsv")
    os.rename("arrayVSfacets.tsv", "main6/arrayVSfacets.tsv")
    os.rename("arrayVSascatNGS.tsv", "main6/arrayVSascatNGS.tsv")
    os.rename("arrayVSsequenza.tsv", "main6/arrayVSsequenza.tsv")


if __name__ == "__main__" :
    main()
