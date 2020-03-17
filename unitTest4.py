#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of introduction of sequenza
"""

"""
The test includes:
    1. Open the sample outputs from sequenza, facets, ascatNGS, and array respectively
    2. Compare all tools against arrays in the whole genome
    3. Get the copy number reported by the 4 tools in specific regions: BRCA1, BRCA2, PALB2, and ATM
"""

import libcomparison as lc
import libgetters as lg
import libstatistics as ls
import libconstants

# Open the output for the sample TCGA-04-1332 output from all the tools
sequenza = lc.convert2region("../9793255c_VS_21fc93b7_Sequenza/TCGA-04-1332_segments.txt", "sequenza")
print(sequenza["purity"])
print(sequenza["ploidy"])
facets = lc.convert2region("../9793255c_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets")
ascat = lc.convert2region("../90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs")
array = lc.convert2region("../73a3a9bb-7dfc-4fc5-9f31-b2630c82010b_Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_F05_437768.grch38.seg.v2.txt", "array")
print("INFO: Arxius oberts satisfactoriament")

print("\nINFO: Resum de les dades obteses en cada eina")
print("Arrays\n-----------\n{}".format(ls.countsXtool(array)))
print("Sequenza\n-----------\n{}".format(ls.countsXtool(sequenza)))
print("FACETS\n-----------\n{}".format(ls.countsXtool(facets)))
print("ascatNGS\n-----------\n{}".format(ls.countsXtool(ascat)))

print("INFO: Comparant les dades amb els arrays")
# TODO mostrar les comparacions entre cadascun dels arrays. Mostrar taules 4x4 i algunes estadistiques
fragments = lc.getFragments(array, sequenza)
tab = lc.doComparison(fragments, array, sequenza)
c1, c2 = ls.calculateCounts(tab)
contingency = ls.doContingency(tab, ["A", "D", "N"])
print(ls.printTable(tab, "Array", "Sequenza", False))


#Regions d'interes. Dades obtingudes des de biogps
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
palb2 = ["16", 23603160, 23641310]
atm = ["11", 108222484, 108369102]

# TODO Consultar summary.md per vore si la mostra te mutacions en BRCA

print("\nINFO: Comparant el LOH en BRCA1")
print("Array output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], array)))
print("Sequenza output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascat)))

print("\nINFO: Comparant el LOH en BRCA2")
print("Array output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], array)))
print("Sequenza output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascat)))

print("\nINFO: Comparant el LOH ATM")
print("Array output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], array)))
print("Sequenza output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], ascat)))

print("\nINFO: Comparant el LOH en PALB2")
print("Array output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], array)))
print("Sequenza output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], ascat)))
