#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of introduction of TCGA ASCAT2
"""

"""
The test includes:
    1. Open example outputs from sequenza, facets, ascatNGS, array, and ASCAT2
    2. Count the aberrations reported by each tool
    3. Get the copy number reported by the all outputs in BRCA1 and BRCA2 genes
"""

import libcomparison as lc
import libgetters as lg
import libstatistics as ls
import main3 as mm

# BRCA1/2 gene coordinates as reported by bioGPS
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]

# Convert the files to REGION format
ascatngs = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs", "error")
sequenza = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_Sequenza/TCGA-04-1332_segments.txt", "sequenza", "error")
facets = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets", "error")
array = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_E11_437726.grch38.seg.v2.txt", "array")
ascat = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/ASCAT2/TCGA-OV.79e63073-7d6d-456b-92c7-a3a7f0216ee7.ascat2.allelic_specific.seg.txt", "ascatarray")

print("INFO: Files opened successfully")
print("INFO: Number of aberrations reported by each tool")
print("Array: {}".format(ls.countsXtool(array)))
print("ASCAT2: {}".format(ls.countsXtool(ascat)))
print("ascatNGS: {}".format(ls.countsXtool(ascatngs)))
print("FACETS: {}".format(ls.countsXtool(facets)))
print("Sequenza: {}".format(ls.countsXtool(sequenza)))
print("----------------------------------------")
print("INFO: Copy number reported in BRCA1")
print("ASCAT2  : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascat)))
print("Array   : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], array)))
print("ascatNGS: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascatngs)))
print("FACETS  : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], facets)))
print("Sequenza: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], sequenza)))
print("----------------------------------------")
print("INFO: Copy number reported in BRCA2")
print("ASCAT2  : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascat)))
print("Array   : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], array)))
print("ascatNGS: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascatngs)))
print("FACETS  : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], facets)))
print("Sequenza: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], sequenza)))
print("----------------------------------------")
print("INFO: Comparing ASCAT vs Array")
regs = lc.getFragments(ascat, array)
comp = lc.doComparison2(regs, ascat, array)
sts = ls.doContingency(comp)
print(ls.printTable(comp, "ASCAT2", "Array", False))
print("\n\tMatthews Correlation Coefficient")
for s in sts :
    print("\t\t{} - {}".format(s, sts[s]["MCC"]))
print("\n\tJaccard Index")
sts = ls.jaccardIndex(comp)
for s in sts :
    print("\t\t{} - {}".format(s, sts[s]))
print("\n\tBase ACC")
print("\t\t{}".format(mm.baseSimilarity(regs, ascat, array)))
print("\n\tRegion ACC")
print("\t\t{}".format(mm.regSimilarity(regs, ascat, array)))
