#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sqlite3

import libcomparison as lc
import libgetters as lg
import libstatistics as ls

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

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
    folder = os.listdir(path)
    for f in folders :
        if f.endswith("_ASCAT") :
            aux = "{abs}/{as}".format(abs = path, as = f)
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

# Buscar els submitters en la base de dades
# Obrir la carpeta d'ASCAT2 i comprovar quants arxius tinc
# Obrir la carpeta d'Arrays i comprovar quants arxius tinc
# Comprovar quantes combinacions tinc per cada eina de LOH

test = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332"
print(getFACETS(test))
print(getAscatNGS(test))
print(getSequenza(test))
"""
print("INFO: Comparing ASCAT2 and Array outputs")
folder1 = "{}/ASCAT2/".format(test)
ascats = os.listdir(folder1)
folder2 = "{}/Array/".format(test)
arrays = os.listdir(folder2)
for a in ascats :
    ascat = lc.convert2region("{}/{}".format(folder1, a), "ascatarray")
    for b in arrays :
        array = lc.convert2region("{}/{}".format(folder2, b), "array")
        regs = lc.getFragments(ascat, array)
        comp = lc.doComparison2(regs, ascat, array)
        sts = ls.doContingency(comp)
        jcc = ls.jaccardIndex(comp)
        num = ls.regionNumber(comp)
        print("{} - {}".format(a[0:8], b[0:8]))
        print("Regions: {}".format(num))
        print("JCC:\tAmp - {}\tDel - {}\tLOH - {}\tNor - {}".format(jcc["A"], jcc["D"], jcc["L"], jcc["N"]))
        print("MCC:\tAmp - {}\tDel - {}\tLOH - {}\tNor - {}".format(sts["A"]["MCC"], sts["D"]["MCC"], sts["L"]["MCC"], sts["N"]["MCC"]))
        print("-------------------------------------------------------------------------")
"""
