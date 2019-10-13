#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Create logR plot for all the samples in the first 81 sample test
"""

"""
Check the compareAscatFACETS folder and get the samples in it
Get the pairs ascatNGS-FACETS in each folder
Get the ploidy and purity from each analysis. Store the data in a txt file
Calculate the
"""

import os
import libextractfile as ext
import libcomparison as comp
import libgetters as ge
import libstatistics as st
import libconstants as cte

path = "/home/ffuster/Desktop/compArrayAscatFacets"
for dirname, dirnames, filenames in os.walk(path) :
    print(dirnames)
    break
