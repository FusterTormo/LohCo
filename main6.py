#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3
import sys

import libcomparison as lc
import libgetters as lg
import libstatistics as ls

"""MAIN PROGRAM
    Compare DNAcopy, ASCAT2, FACETS, ascatNGS, Sequenza, and PURPLE one vs others.
    The output of this main program is the different comparisons done in a tab-delimited text file.
    That means, 36 files: i.e. DNAcopy vs itself, vs ASCAT2, vs FACETS, vs ascatNGS, vs Sequenza, and vs PURPLE.
    Output files are after used in comparaGoldSet and output_vs_5tools R Notebooks to create the plots.
    Stats calculated:
        * Base similarity: Number of bases with same/different aberration reported
        * Region similarity: Number of regions with same/different aberration reported
        * Jaccard Index: Number of bases with the same aberration reported divided by all the bases
        * Matthew's correlation coefficient for each aberration (Amplification, LOH, Deletion and Normal)
        * Jaccard index for each aberration
        * Purity repoted by both tools (if this is the case)
        * Ploidy reported by both tools (if this is the case)
"""

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

def getPurple(path) :
    """Return the PURPLE output files in the submitter folder passed as parameter"""
    purple = []
    folders = os.listdir(path)
    for f in folders :
        if f.endswith("_PURPLE") :
            aux = "{abs}/{fa}/TUMOR.purple.cnv.somatic.tsv".format(abs = path, fa = f)
            if os.path.isfile(aux) :
                purple.append(aux)
    return purple

def getJaccard(comp, a) :
    try :
        return comp[a][a]/(comp[a]["A"] + comp[a]["L"] + comp[a]["N"] + comp[a]["D"] + comp["A"][a] + comp["L"][a] + comp["N"][a] + comp["D"][a] -comp[a][a])
    except ZeroDivisionError:
        return 0


def compareTools(reg1, reg2) :
    regs = lc.getFragments(reg1, reg2)
    comp = lc.doComparison2(regs, reg1, reg2)
    mat = ls.doContingency(comp)
    jcc = ls.jaccardIndex(comp)
    num = ls.regionNumber(regs)
    base = ls.baseSimilarity(regs, reg1, reg2)
    regions = ls.regSimilarity(regs, reg1, reg2)
    jcca = getJaccard(comp, "A")
    jccn = getJaccard(comp, "N")
    jccl = getJaccard(comp, "L")
    jccd = getJaccard(comp, "D")
    st = "{nr}\t{bs}\t{rs}\t{mcca}\t{mccn}\t{mccl}\t{mccd}\t{jcc}\t{jcca}\t{jccn}\t{jccl}\t{jccd}\t{pur1}\t{pur2}\t{plo1}\t{plo2}".format(
        nr = num, bs = base, rs = regions, mcca = mat["A"]["MCC"], mccn = mat["N"]["MCC"], mccl = mat["L"]["MCC"], mccd = mat["D"]["MCC"], jcc = jcc, jcca = jcca, jccn = jccn, jccl = jccl, jccd = jccd, pur1 = reg1["purity"], pur2 = reg2["purity"], plo1 = reg1["ploidy"], plo2 = reg2["ploidy"])
    return st

def createFile(name) :
    """Creates a new file with the name passed as parameter. Prints the header in the file"""
    with open(name, "w") as fi :
        fi.write("submitter\tID1\tID2\tREGIONS\tBASEsim\tREGIONsim\tMCCamp\tMCCnorm\tMCCloh\tMCCdel\tJCC\tJCCamp\tJCCnorm\tJCCloh\tJCCdel\tpurity1\tpurity2\tploidy1\tploidy2\n")

def main() :
    # Constants
    dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
    cancer = "OV"
    cancerpath = "/g/strcombio/fsupek_cancer2/TCGA_bam/"
    if os.path.isdir("main6") :
        print("ERROR: Folder for output already exists. Remove it before to continue")
        sys.exit(1)

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
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(ascat, ascat)))
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
                            fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ascat, arr)))


            # Compare FACETS LOH/CNV outputs with ASCAT2
            facetsFiles = getFACETS(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in facetsFiles :
                    f = lc.convert2region(b, "facets", "error")
                    if not os.path.isfile("ascat2VSfacets.tsv") :
                        createFile("ascat2VSfacets.tsv")
                    with open("ascat2VSfacets.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ascat, f)))

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
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ascat, ngs)))

            # Compare Sequenza LOH/CNV outputs with ASCAT2
            # print("INFO: Comparing ASCAT2 and Sequenza outputs in {}".format(s))
            sequenzaFiles = getSequenza(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in sequenzaFiles :
                    seq = lc.convert2region(b, "sequenza", "error")
                    if not os.path.isfile("ascat2VSsequenza.tsv") :
                        createFile("ascat2VSsequenza.tsv")
                    with open("ascat2VSsequenza.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ascat, seq)))

            purpleFiles = getPurple(workindir)
            for a in ascatFiles :
                ascat = lc.convert2region("{}/{}".format(ascatFolder, a), "ascatarray")
                for b in purpleFiles :
                    purp = lc.convert2region(b, "purple", "error")
                    if not os.path.isfile("ascat2VSpurple.tsv") :
                        createFile("ascat2VSpurple.tsv")
                    with open("ascat2VSpurple.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ascat, purp)))

    # Move the output data to a new folder
    os.mkdir("main6")
    os.rename("ascat2VSascat2.tsv", "main6/ascat2VSascat2.tsv")
    os.rename("ascat2VSarray.tsv", "main6/ascat2VSarray.tsv")
    os.rename("ascat2VSfacets.tsv", "main6/ascat2VSfacets.tsv")
    os.rename("ascat2VSascatNGS.tsv", "main6/ascat2VSascatNGS.tsv")
    os.rename("ascat2VSsequenza.tsv", "main6/ascat2VSsequenza.tsv")
    os.rename("ascat2VSpurple.tsv", "main6/ascat2VSpurple.tsv")

    # Repeat the analysis, but using Arrays as Truth set
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
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                if not os.path.isfile("arrayVSarray.tsv") :
                    createFile("arrayVSarray.tsv")
                with open("arrayVSarray.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(arr, arr)))
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
                            fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(arr, ascat)))


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
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(arr, f)))

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
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(arr, ngs)))

            # Compare Sequenza LOH/CNV outputs with ASCAT2
            # print("INFO: Comparing ASCAT2 and Sequenza outputs in {}".format(s))
            sequenzaFiles = getSequenza(workindir)
            for a in arrayFiles :
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                for b in sequenzaFiles :
                    seq = lc.convert2region(b, "sequenza", "error")
                    if not os.path.isfile("arrayVSsequenza.tsv") :
                        createFile("arrayVSsequenza.tsv")
                    with open("arrayVSsequenza.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(arr, seq)))

            # Compare PURPLE LOH/CNV outputs with DNAcopy
            purpleFiles = getPurple(workindir)
            for a in arrayFiles :
                arr = lc.convert2region("{}/{}".format(arrayFolder, a), "array")
                for b in purpleFiles :
                    purp = lc.convert2region(b, "purple", "error")
                    if not os.path.isfile("arrayVSpurple.tsv") :
                        createFile("arrayVSpurple.tsv")
                    with open("arrayVSpurple.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(arr, purp)))

    os.rename("arrayVSarray.tsv", "main6/arrayVSarray.tsv")
    os.rename("arrayVSascat2.tsv", "main6/arrayVSascat2.tsv")
    os.rename("arrayVSfacets.tsv", "main6/arrayVSfacets.tsv")
    os.rename("arrayVSascatNGS.tsv", "main6/arrayVSascatNGS.tsv")
    os.rename("arrayVSsequenza.tsv", "main6/arrayVSsequenza.tsv")
    os.rename("arrayVSpurple.tsv", "main6/arrayVSpurple.tsv")

    # Repeat the analysis but comparing FACETS vs all the other tools
    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} FACETS".format(s))
        # Get all the FACETS done in the submitter
        facetsFiles = getFACETS(workindir)
        for a in facetsFiles :
            # Compàre FACETS with itself
            f = lc.convert2region(a, "facets", "error")
            if not os.path.isfile("facetsVSfacets.tsv") :
                createFile("facetsVSfacets.tsv")
            with open("facetsVSfacets.tsv", "a") as fi :
                fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(f, f)))

            # Compare with ASCAT2
            ascatFolder = "{}/ASCAT2".format(workindir)
            if os.path.isdir (ascatFolder) :
                # Open ASCAT2 folder and get the files available
                ascatFiles = os.listdir(ascatFolder)
                # Compare ASCAT2 with itself
                for b in ascatFiles :
                    ascat = lc.convert2region("{}/{}".format(ascatFolder, b), "ascatarray", "error")
                    if not os.path.isfile("facetsVSascat2.tsv") :
                        createFile("facetsVSascat2.tsv")
                    with open("facetsVSascat2.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(f, ascat)))

            # Compare with SNP-Arrays
            arrayFolder = "{}/Array".format(workindir)
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                for b in arrayFiles :
                    arr = lc.convert2region("{}/{}".format(arrayFolder, b), "array", "error")
                    if not os.path.isfile("facetsVSarrays.tsv") :
                        createFile("facetsVSarrays.tsv")
                    with open("facetsVSarrays.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(f, arr)))

            # Compare with ascatNGS
            ascatngsFiles = getAscatNGS(workindir)
            for b in ascatngsFiles :
                ngs = lc.convert2region(b, "ascatngs", "error")
                if not os.path.isfile("facetsVSascatNGS.tsv") :
                    createFile("facetsVSascatNGS.tsv")
                with open("facetsVSascatNGS.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(f, ngs)))

            # Compare with Sequenza
            sequenzaFiles = getSequenza(workindir)
            for b in sequenzaFiles :
                seq = lc.convert2region(b, "sequenza", "error")
                if not os.path.isfile("facetsVSsequenza.tsv") :
                    createFile("facetsVSsequenza.tsv")
                with open("facetsVSsequenza.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(f, seq)))

            # Compare with PURPLE
            purpleFiles = getPurple(workindir)
            for b in purpleFiles :
                purp = lc.convert2region(b, "purple", "error")
                if not os.path.isfile("facetsVSpurple.tsv") :
                    createFile("facetsVSpurple.tsv")
                with open("facetsVSpurple.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(f, purp)))

    os.rename("facetsVSfacets.tsv", "main6/facetsVSfacets.tsv")
    os.rename("facetsVSascat2.tsv", "main6/facetsVSascat2.tsv")
    os.rename("facetsVSarrays.tsv", "main6/facetsVSarrays.tsv")
    os.rename("facetsVSascatNGS.tsv", "main6/facetsVSascatNGS.tsv")
    os.rename("facetsVSsequenza.tsv", "main6/facetsVSsequenza.tsv")
    os.rename("facetsVSpurple.tsv", "main6/facetsVSpurple.tsv")

    # Repeat the analysis, but comparing ascatNGS vs all the other tools
    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} ascatNGS".format(s))
        # Get all the ascatNGS done in the submitter
        ascatngsFiles = getAscatNGS(workindir)
        for a in ascatngsFiles :
            # Compare ascatNGS vs itself
            ngs = lc.convert2region(a, "ascatngs", "error")
            if not os.path.isfile("ascatNGSVSascatNGS.tsv") :
                createFile("ascatNGSVSascatNGS.tsv")
            with open("ascatNGSVSascatNGS.tsv", "a") as fi :
                fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(ngs, ngs)))

            # Compare with ASCAT2
            ascatFolder = "{}/ASCAT2/".format(workindir)
            if os.path.isdir(ascatFolder) :
                ascatFiles = os.listdir(ascatFolder)
                for b in ascatFiles :
                    ascat = lc.convert2region("{}{}".format(ascatFolder, b), "ascatarray", "error")
                    if not os.path.isfile("ascatNGSVSascat2.tsv") :
                        createFile("ascatNGSVSascat2.tsv")
                    with open("ascatNGSVSascat2.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ngs, ascat)))

            # Compare with SNP-Arrays
            arrayFolder = "{}/Array/".format(workindir)
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                for b in arrayFiles :
                    arr = lc.convert2region("{}{}".format(arrayFolder, b), "array", "error")
                    if not os.path.isfile("ascatNGSVSarrays.tsv") :
                        createFile("ascatNGSVSarrays.tsv")
                    with open("ascatNGSVSarrays.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ngs, arr)))

            # Compare with FACETS
            facetsFiles = getFACETS(workindir)
            for b in facetsFiles :
                f = lc.convert2region(b, "facets", "error")
                if not os.path.isfile("ascatNGSVSfacets.tsv") :
                    createFile("ascatNGSVSfacets.tsv")
                with open("ascatNGSVSfacets.tsv", "a") as fi :

                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ngs, f)))

            # Compare with Sequenza
            sequenzaFiles = getSequenza(workindir)
            for b in sequenzaFiles :
                seq = lc.convert2region(b, "sequenza", "error")
                if not os.path.isfile("ascatNGSVSsequenza.tsv") :
                    createFile("ascatNGSVSsequenza.tsv")
                with open("ascatNGSVSsequenza.tsv", "a") as fi:
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ngs, seq)))

            # Compare with PURPLE
            purpleFiles = getPurple(workindir)
            for b in purpleFiles :
                purp = lc.convert2region(b, "purple", "error")
                if not os.path.isfile("ascatNGSVSpurple.tsv") :
                    createFile("ascatNGSVSpurple.tsv")
                with open("ascatNGSVSpurple.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(ngs, purp)))

    os.rename("ascatNGSVSascat2.tsv", "main6/ascatNGSVSascat2.tsv")
    os.rename("ascatNGSVSarrays.tsv", "main6/ascatNGSVSarrays.tsv")
    os.rename("ascatNGSVSfacets.tsv", "main6/ascatNGSVSfacets.tsv")
    os.rename("ascatNGSVSascatNGS.tsv", "main6/ascatNGSVSascatNGS.tsv")
    os.rename("ascatNGSVSsequenza.tsv", "main6/ascatNGSVSsequenza.tsv")
    os.rename("ascatNGSVSpurple.tsv", "main6/ascatNGSVSpurple.tsv")

    # Repeat the analysis, but comparing Sequenza vs all the other approximations
    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} Sequenza".format(s))
        # Get all the Sequenza done in the submitter
        sequenzaFiles = getSequenza(workindir)
        for a in sequenzaFiles :
            # Compare sequenza vs itself
            seq = lc.convert2region(a, "sequenza", "error")
            if not os.path.isfile("sequenzaVSsequenza.tsv") :
                createFile("sequenzaVSsequenza.tsv")
            with open("sequenzaVSsequenza.tsv", "a") as fi :
                fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(seq, seq)))

            # Compare with ASCAT2
            ascatFolder = "{}/ASCAT2".format(workindir)
            if os.path.isdir(ascatFolder) :
                ascatFiles = os.listdir(ascatFolder)
                for b in ascatFiles :
                    ascat = lc.convert2region("{}/{}".format(ascatFolder, b), "ascatarray", "error")
                    if not os.path.isfile("sequenzaVSascat2.tsv") :
                        createFile("sequenzaVSascat2.tsv")
                    with open("sequenzaVSascat2.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(seq, ascat)))

            # Compare with SNP-Arrays
            arrayFolder = "{}/Array".format(workindir)
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                for b in arrayFiles :
                    arr = lc.convert2region("{}/{}".format(arrayFolder, b), "array", "error")
                    if not os.path.isfile("sequenzaVSarrays.tsv") :
                        createFile("sequenzaVSarrays.tsv")
                    with open("sequenzaVSarrays.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(seq, arr)))

            # Compare with FACETS
            facetsFiles = getFACETS(workindir)
            for b in facetsFiles :
                f = lc.convert2region(b, "facets", "error")
                if not os.path.isfile("sequenzaVSfacets.tsv") :
                    createFile("sequenzaVSfacets.tsv")
                with open("sequenzaVSfacets.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(seq, f)))

            # Compare with ascatNGS
            ascatngsFiles = getAscatNGS(workindir)
            for b in ascatngsFiles :
                ngs = lc.convert2region(b, "ascatngs", "error")
                if not os.path.isfile("sequenzaVSascatNGS.tsv") :
                    createFile("sequenzaVSascatNGS.tsv")
                with open("sequenzaVSascatNGS.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(seq, ngs)))

            # Compare with PURPLE
            purpleFiles = getPurple(workindir)
            for b in purpleFiles :
                purp = lc.convert2region(b, "purple", "error")
                if not os.path.isfile("sequenzaVSpurple.tsv") :
                    createFile("sequenzaVSpurple.tsv")
                with open("sequenzaVSpurple.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(seq, purp)))

    os.rename("sequenzaVSascat2.tsv", "main6/sequenzaVSascat2.tsv")
    os.rename("sequenzaVSarrays.tsv", "main6/sequenzaVSarrays.tsv")
    os.rename("sequenzaVSfacets.tsv", "main6/sequenzaVSfacets.tsv")
    os.rename("sequenzaVSascatNGS.tsv", "main6/sequenzaVSascatNGS.tsv")
    os.rename("sequenzaVSsequenza.tsv", "main6/sequenzaVSsequenza.tsv")
    os.rename("sequenzaVSpurple.tsv", "main6/sequenzaVSpurple.tsv")

    # Repeat the analysis but comparing PURPLE vs all the other tools
    for sub in submitters :
        s = sub[0]
        workindir = "{}/{}/{}".format(cancerpath, cancer, s)
        print("INFO: Checking {} PURPLE".format(s))
        # Get all the Sequenza done in the submitter
        purpleFiles = getPurple(workindir)
        for a in purpleFiles :
            # Compare PURPLE vs itself
            purp = lc.convert2region(a, "purple", "error")
            if not os.path.isfile("purpleVSpurple.tsv") :
                createFile("purpleVSpurple.tsv")
            with open("purpleVSpurple.tsv", "a") as fi :
                fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = a, cmp = compareTools(purp, purp)))

            # Compare with ASCAT2
            ascatFolder = "{}/ASCAT2".format(workindir)
            if os.path.isdir(ascatFolder) :
                ascatFiles = os.listdir(ascatFolder)
                for b in ascatFiles :
                    ascat = lc.convert2region("{}/{}".format(ascatFolder, b), "ascatarray", "error")
                    if not os.path.isfile("purpleVSascat2.tsv") :
                        createFile("purpleVSascat2.tsv")
                    with open("purpleVSascat2.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(purp, ascat)))

            # Compare with SNP-Arrays
            arrayFolder = "{}/Array".format(workindir)
            if os.path.isdir(arrayFolder) :
                arrayFiles = os.listdir(arrayFolder)
                for b in arrayFiles :
                    arr = lc.convert2region("{}/{}".format(arrayFolder, b), "array", "error")
                    if not os.path.isfile("purpleVSarrays.tsv") :
                        createFile("purpleVSarrays.tsv")
                    with open("purpleVSarrays.tsv", "a") as fi :
                        fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(purp, arr)))

            # Compare with FACETS
            facetsFiles = getFACETS(workindir)
            for b in facetsFiles :
                f = lc.convert2region(b, "facets", "error")
                if not os.path.isfile("purpleVSfacets.tsv") :
                    createFile("purpleVSfacets.tsv")
                with open("purpleVSfacets.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(purp, f)))

            # Compare with ascatNGS
            ascatngsFiles = getAscatNGS(workindir)
            for b in ascatngsFiles :
                ngs = lc.convert2region(b, "ascatngs", "error")
                if not os.path.isfile("purpleVSascatNGS.tsv") :
                    createFile("purpleVSascatNGS.tsv")
                with open("purpleVSascatNGS.tsv", "a") as fi :
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(purp, ngs)))

            # Compare with Sequenza
            sequenzaFiles = getSequenza(workindir)
            for b in sequenzaFiles :
                seq = lc.convert2region(b, "sequenza", "error")
                if not os.path.isfile("purpleVSsequenza.tsv") :
                    createFile("purpleVSsequenza.tsv")
                with open("purpleVSsequenza.tsv", "a") as fi:
                    fi.write("{sub}\t{id1}\t{id2}\t{cmp}\n".format(sub = sub[0], id1 = a, id2 = b, cmp = compareTools(purp, seq)))

    os.rename("purpleVSascat2.tsv", "main6/purpleVSascat2.tsv")
    os.rename("purpleVSarrays.tsv", "main6/purpleVSarrays.tsv")
    os.rename("purpleVSfacets.tsv", "main6/purpleVSfacets.tsv")
    os.rename("purpleVSascatNGS.tsv", "main6/purpleVSascatNGS.tsv")
    os.rename("purpleVSsequenza.tsv", "main6/purpleVSsequenza.tsv")
    os.rename("purpleVSpurple.tsv", "main6/purpleVSpurple.tsv")

if __name__ == "__main__" :
    main()
