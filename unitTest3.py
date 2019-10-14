#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Create logR plot for all the samples in the first 81 sample test
"""

"""
Check the compArrayAscatFACETS folder and get the samples in it
Get the arrays in each sample
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

def launchAnalysis(folder, array, ascat, facets) :
    pythonpath = os.path.dirname(os.path.realpath(__file__))
    #Go to the working directory to run the comparison analysis
    print("INFO: Analysing {}".format(folder))
    os.chdir(folder)
    ascatReg = None
    facetsReg = None
    arrayReg = None
    arrayPath = "../{}".format(array)
    id = folder.split("/")[-1]
    if ascat != "" :
        ascatReg = comp.convert2region(ascat, "ascat")
    if facets != "" :
        facetsReg = comp.convert2region(facets, "facets")

    arrayReg = comp.convert2region(arrayPath, "array")

    if ascatReg != None and facetsReg != None : #Calculate logR between AscatNGS and FACETS
        regA_F = comp.getFragments(ascatReg, facetsReg)
        st.logRcomp(regA_F, ascatReg, facetsReg, "ascatNGS", "FACETS")

    if facetsReg != None and arrayReg != None :
        regAr_F = comp.getFragments(arrayReg, facetsReg)
        mt1 = comp.doComparison2(regAr_F, facetsReg, arrayReg)
        jc1 = st.jaccardIndex(mt1, ["A", "D"])
        mt1 = comp.doComparison(regAr_F, facetsReg, arrayReg)
        cm1 = st.doContingency(mt1, ["A", "D"])

    if ascatReg != None and arrayReg != None :
        regAr_A = comp.getFragments(arrayReg, ascatReg)
        mt2 = comp.doComparison2(regAr_F, ascatReg, arrayReg)
        jc2 = st.jaccardIndex(mt2, ["A", "D"])
        mt2 = comp.doComparison(regAr_F, ascatReg, arrayReg)
        cm2 = st.doContingency(mt2, ["A", "D"])

    #Return to current python path
    os.chdir(pythonpath)
    #Store the summary data in the corresponding files
    if ascatReg != None :
        with open("ascatPurities.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, ascatReg["purity"]))
        with open("ascatPloidies.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, ascatReg["ploidy"]))
        with open("ascatJaccardAmplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, jc2["A"]))
        with open("ascatJaccardDeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, jc2["D"]))
        with open("ascatACCamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["A"]["ACC"]))
        with open("ascatTPRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["A"]["TPR"]))
        with open("ascatTNRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["A"]["TNR"]))
        with open("ascatPPVamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["A"]["PPV"]))
        with open("ascatFDRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["A"]["ACC"]))
        with open("ascatACCdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["D"]["ACC"]))
        with open("ascatTPRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["D"]["TPR"]))
        with open("ascatTNRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["D"]["TNR"]))
        with open("ascatPPVdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["D"]["PPV"]))
        with open("ascatFDRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm2["D"]["ACC"]))
    if facetsReg != None :
        with open("facetsPurities.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, facetsReg["purity"]))
        with open("facetsPloidies.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, facetsReg["ploidy"]))
        with open("facetsJaccardAmplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, jc1["A"]))
        with open("facetsJaccardDeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, jc1["D"]))
        with open("facetsACCamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["A"]["ACC"]))
        with open("facetsTPRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["A"]["TPR"]))
        with open("facetsTNRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["A"]["TNR"]))
        with open("facetsPPVamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["A"]["PPV"]))
        with open("facetsFDRamplification.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["A"]["ACC"]))
        with open("facetsACCdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["D"]["ACC"]))
        with open("facetsTPRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["D"]["TPR"]))
        with open("facetsTNRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["D"]["TNR"]))
        with open("facetsPPVdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["D"]["PPV"]))
        with open("facetsFDRdeletion.txt", "a") as fi :
            fi.write("{}\t{}\n".format(id, cm1["D"]["ACC"]))

#Get the samples in the directory. Each sample has its own directory
path = "/home/ffuster/Desktop/compArrayAscatFacets"
facetsStruct = "facets_comp_cncf.tsv"
ascatStruct = "copynumber.caveman.csv"
for dirname, dirnames, filenames in os.walk(path) :
    dirs = dirnames
    break

print("INFO: Analysis will be done in {} samples".format(len(dirs)))
#For each directory,
for d in dirs :
    current = "{}/{}".format(path, d)
    for dirname, dirnames, filenames in os.walk(current) :
        analyses = dirnames #Get the analyses tumor_vs_control done in the sample
        arrays = filenames #Get the arrays that each sample has
        for a in analyses :
            subpath = "{}/{}".format(current, a)
            for aux, aux2, fils in os.walk(subpath) :
                break
            #At this step we have the array, ascat and facets done in folder 'a'. Start the analyses
            ascat = ""
            facets = ""
            for f in fils :
                if f == facetsStruct :
                    facets = f
                if f.endswith(ascatStruct) :
                    ascat = f
            launchAnalysis(subpath, arrays[0], ascat, facets) # TODO hi ha una mostra que te 3 arrays
        break
