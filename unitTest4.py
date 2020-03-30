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

# Open the output for the sample TCGA-04-1332 output from all the tools
sequenza = lc.convert2region("../9793255c_VS_f4b549d0_Sequenza/TCGA-04-1332_segments.txt", "sequenza")
facets = lc.convert2region("../9793255c_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets")
ascat = lc.convert2region("../90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs")
array = lc.convert2region("../73a3a9bb-7dfc-4fc5-9f31-b2630c82010b_Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_F05_437768.grch38.seg.v2.txt", "array")
print("INFO: Arxius oberts satisfactoriament")

# Print the counts in each file
print("\nINFO: Resum de les dades obteses en cada eina")
car = ls.countsXtool(array)
cs = ls.countsXtool(sequenza)
cf = ls.countsXtool(facets)
cas = ls.countsXtool(ascat)
print("\t A\t| N \t| L \t| D")
print("Array \t {amp}\t| {normal}\t| {loh}\t| {dele}".format(amp = car["A"], normal = car["N"], loh = car["L"], dele = car["D"]))
print("FACETS\t {amp}\t| {normal}\t| {loh}\t| {dele}".format(amp = cf["A"], normal = cf["N"], loh = cf["L"], dele = cf["D"]))
print("ascatN\t {amp}\t| {normal}\t| {loh}\t| {dele}".format(amp = cas["A"], normal = cas["N"], loh = cas["L"], dele = cas["D"]))
print("Sequen\t {amp}\t| {normal}\t| {loh}\t| {dele}".format(amp = cs["A"], normal = cs["N"], loh = cs["L"], dele = cs["D"]))

print("\nINFO: Comparant les dades amb els arrays")
fragments = lc.getFragments(array, facets)
tab = lc.doComparison(fragments, array, facets)
c1, c2 = ls.calculateCounts(tab)
cAvF = ls.doContingency(tab, ["A", "D", "N"])
tabAvF = ls.printTable(tab, "Array", "FACETS", False).split("\n")

fragments = lc.getFragments(array, ascat)
tab = lc.doComparison(fragments, array, ascat)
cAvA = ls.doContingency(tab, ["A", "D", "N"])
tabAvA = ls.printTable(tab, "Array", "ascatNGS", False).split("\n")

fragments = lc.getFragments(array, sequenza)
tab = lc.doComparison(fragments, array, sequenza)
cAvS = ls.doContingency(tab, ["A", "D", "N"])
tabAvS = ls.printTable(tab, "Array", "Sequenza", False).split("\n")
print("{}\t\t\t\t||{}\t\t\t||{}".format(tabAvF[0], tabAvA[0], tabAvS[0]))
print("{}\t\t||{}\t\t||{}".format(tabAvF[1], tabAvA[1], tabAvS[1]))
print("{}\t\t||{}\t\t||{}".format(tabAvF[2], tabAvA[2], tabAvS[2]))
print("{}\t\t||{}\t\t||{}".format(tabAvF[3], tabAvA[3], tabAvS[3]))
print("{}\t\t||{}\t\t||{}".format(tabAvF[4], tabAvA[4], tabAvS[4]))
print("{}\t\t||{}\t\t||{}".format(tabAvF[5], tabAvA[5], tabAvS[5]))
print("=====================================================================")
print("TODO: Add information regarding the contingency table for the alterations (A, D, N)")


# Get the interesting information from the variant calling summary file
summaryFile = "../summary.md"
data = {}
with open(summaryFile, "r") as fi :
	for l in fi :
        if l.startswith("TCGA-04-1332") :
            aux = l.split("|")
            if (aux[1].startswith("9793") or aux[1].startswith("21fc")) and aux[3] == "Platypus" :
                aux2 = aux[6].strip().split("&rarr;")
                data[aux[1]] = {"variant" : aux2[1], "gene" : aux2[0]}

#Regions of interest. Data extracted from biogps
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
palb2 = ["16", 23603160, 23641310]
atm = ["11", 108222484, 108369102]

print("|Cas|Variant|     |Tipus|FACETS|    |   |     |ascatNGS|  |   |     |Sequenza|  |   |     |")
print("|Cas|Tumor|Control|Tipus|BRCA1|BRCA2|ATM|PALB2|BRCA1|BRCA2|ATM|PALB2|BRCA1|BRCA2|ATM|PALB2|")
print("|---|-----|-------|-----|-----|-----|---|-----|-----|-----|---|-----|-----|-----|---|-----|")
print("TCGA-04-1332|{tm}|{cn}|{type}|{f1}|{f2}|{f3}|{f4}|{a1}|{a2}|{a3}|{a4}|{s1}|{s2}|{s3}|{s4}|".format(tm = , cn = , type = "?",
    f1 = lg.getCopyNumber(brca1[1:3], brca1[0], facets), f2 = lg.getCopyNumber(brca2[1:3], brca2[0], facets), f3 = lg.getCopyNumber(atm[1:3], atm[0], facets),
    f4 = lg.getCopyNumber(palb2[1:3], palb2[0], facets),
    a1 = lg.getCopyNumber(brca1[1:3], brca1[0], ascat), a2 = lg.getCopyNumber(brca2[1:3], brca2[0], ascat), a3 = lg.getCopyNumber(atm[1:3], atm[0], ascat),
    a4 = lg.getCopyNumber(palb2[1:3], palb2[0], ascat),
    s1 = lg.getCopyNumber(brca1[1:3], brca1[0], sequenza), s2 = lg.getCopyNumber(brca2[1:3], brca2[0], sequenza), s3 = lg.getCopyNumber(atm[1:3], atm[0], sequenza),
    s4 = lg.getCopyNumber(palb2[1:3], palb2[0], sequenza)))
# print("\nINFO: Comparant el LOH en BRCA1")
# print("Array output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], array)))
# print("Sequenza output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], sequenza)))
# print("FACETS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], facets)))
# print("ascatNGS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascat)))
#
# print("\nINFO: Comparant el LOH en BRCA2")
# print("Array output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], array)))
# print("Sequenza output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], sequenza)))
# print("FACETS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], facets)))
# print("ascatNGS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascat)))
#
# print("\nINFO: Comparant el LOH ATM")
# print("Array output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], array)))
# print("Sequenza output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], sequenza)))
# print("FACETS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], facets)))
# print("ascatNGS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], ascat)))
#
# print("\nINFO: Comparant el LOH en PALB2")
# print("Array output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], array)))
# print("Sequenza output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], sequenza)))
# print("FACETS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], facets)))
# print("ascatNGS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], ascat)))
