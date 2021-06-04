#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Testing the function that counts the number of bases that have each aberration
"""

import libcomparison as lc
import libstatistics as ls

print("INFO: Loading example from FACETS")
f = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/")
print("INFO: Loading example from ascatNGS")
n = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/")
print("INFO: Loading example from Sequenza")
s = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/")
print("INFO: Loading example from PURPLE")
p = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/")
print("INFO: Number of regions that have reported each tool")
print(ls.countsXtool(f))
print(ls.countsXtool(a))
print(ls.countsXtool(s))
print(ls.countsXtool(p))
print("INFO: Number of bases that have reported each tool")
print(ls.basesXtool(f))
print(ls.basesXtool(a))
print(ls.basesXtool(s))
print(ls.basesXtool(p))
