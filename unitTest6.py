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
f = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-13-0757/1f1f7441_VS_26fa0e90_FACETS/facets_comp_cncf.tsv", "facets")
print("INFO: Loading example from ascatNGS")
a = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-13-0757/1f1f7441_VS_26fa0e90_ASCAT/TCGA-04-1331-01A-01W.copynumber.caveman.csv", "ascatngs")
print("INFO: Loading example from Sequenza")
s = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-13-0757/1f1f7441_VS_26fa0e90_Sequenza/TCGA-04-1331_segments.txt", "sequenza")
print("INFO: Loading example from PURPLE")
p = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-13-0757/1f1f7441_VS_26fa0e90_PURPLE/TUMOR.purple.cnv.somatic.tsv", "purple")
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
fullLength = sum(ct.lenChromosomes.values())
print(ct.lenChromosomes.values())
print("INFO: Number of bases that have reported each tool")
ff = ls.basesXtool(f)
print("FACETS")
print(ff)
print("Mean CN: {}\n".format(ls.meanCoverage(f)))
print("ASCAT2")
aa = ls.basesXtool(a)
print(aa)
print("Mean CN: {}\n".format(ls.meanCoverage(a)))
print("Sequenza")
ss = ls.basesXtool(s)
print(ss)
print("Mean CN: {}\n".format(ls.meanCoverage(s)))
print("PURPLE")
pp = ls.basesXtool(p)
print("Mean CN: {}".format(ls.meanCoverage(p)))
