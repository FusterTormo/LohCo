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
    Estadistiques?????
"""

import libcomparison as lc
import libgetters as lg

sequenza = lc.convert2region("../9793255c_VS_21fc93b7_Sequenza/TCGA-04-1332_segments.txt", "sequenza")
facets = lc.convert2region("../9793255c_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets")
ascat = lc.convert2region("../90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs")
array = lc.convert2region("", "array")
print("INFO: Arxius oberts satisfactoriament")
print(array)

print("INFO: Comparant les dades amb els arrays")

# TODO mostrar les comparacions entre cadascun dels arrays. Mostrar taules 4x4 i algunes estadistiques
# IDEA: usar lc.getFragments primer

#Regions d'interes. Dades obtingudes des de biogps
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
palb2 = ["16", 23603160, 23641310]
atm = ["11", 108222484, 108369102]

# TODO Consultar summary.md per vore si la mostra te mutacions en BRCA

print("\nINFO: Comparant el LOH en BRCA1")
print("Sequenza output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascat)))

print("\nINFO: Comparant el LOH en BRCA2")
print("Sequenza output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascat)))

print("\nINFO: Comparant el LOH ATM")
print("Sequenza output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(atm[1:3], atm[0], ascat)))

print("\nINFO: Comparant el LOH en PALB2")
print("Sequenza output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], sequenza)))
print("FACETS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], facets)))
print("ascatNGS output: {}".format(lg.getCopyNumber(palb2[1:3], palb2[0], ascat)))