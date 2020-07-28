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

# Buscar els submitters en la base de dades
# Obrir la carpeta d'ASCAT2 i comprovar quants arxius tinc
# Obrir la carpeta d'Arrays i comprovar quants arxius tinc
# Comprovar quantes combinacions tinc per cada eina de LOH

test = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332"
folder1 = "{}/ASCAT2/".format(test)
ascats = os.listdir(folder1)
folder2 = "{}/Array/".format(test)
arrays = os.listdir(folder2)
print(ascats)
print(arrays)
for a in ascats :
    ascat = lc.convert2region("{}/{}".format(folder1, a))
    for b in arrays :
        array = lc.convert2region("{}/{}".format(folder2, b))
        regs = lc.getFragments(ascat, array)
        comp = lc.doComparison2(regs, ascat, array)
        sts = ls.doContingency(comp)
        print(sts["A"])
