#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of comparison between a FACETS example vs ascatNGS example
"""

"""
The test includes :
    1) Read FACETS and ascatNGS inputs
    2) Fragment the regions to get same coordinates in each
    3) Get ploidy and purity from each output
    4) Plot the correlation between logRs
    5) ggplot of total copy number (tcn)
    6) ggplot of minor copy number (lcn)
"""

import libextractfile as exfi
import libcomparison as compi
import libgetters as ge
import libstatistics as sts
import libconstants as cts
import os

print "INFO: Test unitario para comprobar graficas comparando un ejemplo de ascatNGS y uno de FACETS"
print "TEST 1) Extraer los datos"
ascat = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/TCGA-13-0887-01A-01W.copynumber.caveman.csv", "ascat")
facets = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/facets_comp_cncf.tsv", "FACETS")
print "TEST 2) Dividir las regiones para obtener regiones en comun"
regs = compi.getFragments(facets, ascat)
print "TEST 3) Dibujar la concordancia entre los logR"
sts.logRcomp(regs, facets, ascat, "FACETS", "ASCAT")
print "TEST 4) Dibujar los copy number counts usando la libreria ggplot"
sts.doGGplotFiles(facets, ascat, "FACETS", "ASCAT")
print "INFO: Limpiando..."
os.remove("FACETS_ASCAT_lcn.png")
os.remove("FACETS_ASCAT_logRcomp.tsv")
os.remove("FACETS_ASCAT_tcn.png")
os.remove("lcn_ggplot.tsv")
os.remove("logRcorrelation.png")
os.remove("Rplots.pdf")
os.remove("TCN_ggplot.tsv")
