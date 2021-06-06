#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Testing the function that counts the number of bases that have each aberration
"""

import libcomparison as lc
import libstatistics as ls
import libconstants as ct
import libgetters as lg

print("INFO: Loading example from FACETS")
f = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1331/c21ab280_VS_82704a7d_FACETS/facets_comp_cncf.tsv", "facets")
print("INFO: Loading example from ascatNGS")
a = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1331/c21ab280_VS_82704a7d_ASCAT/TCGA-04-1331-01A-01W.copynumber.caveman.csv", "ascatngs")
print("INFO: Loading example from Sequenza")
s = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1331/c21ab280_VS_82704a7d_Sequenza/TCGA-04-1331_segments.txt", "sequenza")
print("INFO: Loading example from PURPLE")
p = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1331/c21ab280_VS_82704a7d_PURPLE/TUMOR.purple.cnv.somatic.tsv", "purple")
allbases = []
it = 0
current = 0



"""
print("INFO: Number of regions that have reported each tool")
print(ls.countsXtool(f))
print(ls.countsXtool(a))
print(ls.countsXtool(s))
print(ls.countsXtool(p))
"""
print("INFO: Number of bases that have reported each tool")
print(ls.basesXtool(f))
print(ls.basesXtool(a))
print(ls.basesXtool(s))
print(ls.basesXtool(p))
