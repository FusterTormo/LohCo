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

print("INFO: Test unitario para comprobar graficas comparando un ejemplo de ascatNGS y uno de FACETS")
print("TEST 1) Extraer los datos")
ascat = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/TCGA-09-0369/TCGA-09-0369_40e311a4_VS_f4441d6e/H_GP-09-0369-01A-01W-0372-09-1.copynumber.caveman.csv", "ascat")
facets = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/TCGA-09-0369/TCGA-09-0369_40e311a4_VS_f4441d6e/facets_comp_cncf.tsv", "FACETS")
print("TEST 2) Dividir las regiones para obtener regiones en comun")
regs = compi.getFragments(facets, ascat)
print("TEST 3) Dibujar la concordancia entre los logR")
try :
    sts.logRcomp(regs, facets, ascat, "FACETS", "ASCAT")
except ValueError :
    print("ERROR: Cannot create the logR plot")
print("TEST 4) Dibujar los copy number counts usando la libreria ggplot")
sts.doGGplotFiles(facets, ascat, "FACETS", "ASCAT")
print("TEST 5) Crear un bed con las regiones reportadas por cada archivo y las regiones en comun")
sts.print2Bed(ascat, "ascatNGS", facets, "FACETS", regs)
print("INFO: Limpiando...")
if os.path.isfile("FACETS_ASCAT_lcn.png") :
    os.remove("FACETS_ASCAT_lcn.png")
if os.path.isfile("FACETS_ASCAT_logRcomp.tsv") :
    os.remove("FACETS_ASCAT_logRcomp.tsv")
if os.path.isfile("FACETS_ASCAT_tcn.png") :
    os.remove("FACETS_ASCAT_tcn.png")
if os.path.isfile("lcn_ggplot.tsv") :
    os.remove("lcn_ggplot.tsv")
if os.path.isfile("logRcorrelation.png") :
    os.remove("logRcorrelation.png")
if os.path.isfile("Rplots.pdf") :
    os.remove("Rplots.pdf")
if os.path.isfile("TCN_ggplot.tsv") :
    os.remove("TCN_ggplot.tsv")
if os.path.isfile("ascatNGS_FACETS_all_regions.bed") :
    os.remove("ascatNGS_FACETS_all_regions.bed")
