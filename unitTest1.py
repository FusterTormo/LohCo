#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of comparison between a FACETS example against the "silver standard" CNV array obtained from TCGA
"""

"""
The test includes:
    1) Read output from array example, and from FACETS example
    2) Fragment the regions to get same coordinates in each
    3) Calculate the 4x4 table
    4) Extract confusion matrix for (A)mplification, (D)eletion, and (N)ormal copy number
    5) Calculate the Jaccard index for the same aberrations
"""

import libextractfile as exfi
import libcomparison as compi
import libgetters as ge
import libstatistics as sts
import libconstants as cts

print "INFO: Test unitario para comparar el output de un ejemplo de FACETS con los datos del array descargado desde TCGA"
print "TEST 1) Extraer datos"
ar = compi.convert2region("input_examples/arrays/TCGA-13-0887_tumor.txt", "array")
fa = compi.convert2region("input_examples/facets_comp_cncf.tsv", "FACETS")
print "TEST 2) Buscando las regiones en comun para estudio"
regs = compi.getFragments(ar, fa)
print "TEST 3) Crear la tabla comparativa 4x4"
dc = compi.doComparison(regs, ar, fa)
print sts.printTable(dc, "Array", "FACETS", False)
print "TEST 4) Resultados de la matriz de confusion para cada aberracion"
c1, c2 = sts.calculateCounts(dc)
dicContingency = sts.doContingency(dc, c1, c2, ["A", "D", "N"])
print "\tAmplificacion\n\t{}\n\n".format(dicContingency["A"])
print "\tDelecion\n\t{}\n\n".format(dicContingency["D"])
print "\tNormal\n\t{}\n".format(dicContingency["N"])
print "TEST 5) Jaccard index de cada una de las aberraciones"
jci2 = compi.doComparison2(regs, ar, fa)
dicJaccard = sts.jaccardIndex(jci2, ["A", "D", "N"])
print "\tAmplificacion\n\t{}\n\n".format(dicJaccard["A"])
print "\tDelecion\n\t{}\n\n".format(dicJaccard["D"])
print "\tNormal\n\t{}\n".format(dicJaccard["N"])
