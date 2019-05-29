#!/usr/bin/python
# -*- coding: utf-8 -*-
import libcomparison as lc
import sys

def extractFromFile(tip, path) :
    tip = tip.lower()
    if tip == "array" :
        (reg, cn) = lc.extractArray(path)
    elif tip == "facets" :
        (reg, cn) = lc.extractFacets(path)
    elif tip == "ascat" or tip == "ascatngs" :
        (reg, cn) = lc.extractAscat(path)
    elif tip == "sequenza" :
        print "ERROR: Function pending to finish"
        sys.exit()
    elif tip == "titan" :
        print "ERROR: Function pending to finish"
        sys.exit()
    elif tip == "purple" :
        print "ERROR: Function pending to finish"
        sys.exit()
    else :
        raise IOError("ERROR: Type of file not found. Cannot continue")

    return (reg, cn)

def main() :
    if len(sys.argv) != 5 :
        #TODO continuar amb l'ajuda en cas de parametres escrits malament
        print "ERROR: Number of parameters is not correct"
        print "USAGE: compareLOH.py path_to_output_file1 type_of_file1 path_to_output_file2 type_of_file2"
        print "Where:\n\tpath_to_file is the path to the output file[...]"
        print "\ttype_of_file can be FACETS, ascatNGS, Array [...]"
    else :
        sw1 = sys.argv[2]
        sw2 = sys.argv[4]
        (reg1, cn1) = extractFromFile(sw1, sys.argv[1])
        (reg2, cn2) = extractFromFile(sw2, sys.argv[3])

        fgm = lc.getFragments(reg1, reg2)
        lc.print2Bed(cn1, sw1, cn2, sw2, fgm)
        lc.checkCopyNumber(fgm, cn1, cn2, sw1[0:2].upper(), sw2[0:2].upper())

    """
    OLD main program
    (ar, cn_ar) = lc.extractArray("arrays/dummy.txt")
    (fa, cn_fa) = lc.extractFacets("facets_chr20.tsv")
    print "ARRAY calls\n{}\n".format(ar)
    print "FACETS calls\n{}\n".format(fa)
    #TODO review if there is a value that is interpreted as amplification/deletion/normal and the value of tcn is not in accordance
    #print cn_fa

    #(ar, cn) = lc.extractArray("arrays/TCGA-13-0887_tumor.txt")
    #(fa, cfa) = lc.extractFacets("facets_comp_cncf.tsv")

    fgm = lc.getFragments(ar, fa)
    print "All fragments\n{}\n".format(fgm)

    lc.checkCopyNumber(fgm, cn_ar, cn_fa)
    """

main()
