#!/usr/bin/python
# -*- coding: utf-8 -*-

import libconstants as lc

"""
MAIN: All functions to return data from a REGION variable
"""

"""
FUNCTIONS:
    - getRegion : Return the start and end position from the chr list
    - getCopyNumber : Return the copy number classification for the REGION passed as parameter in the region passed as parameter
    - getTotalCN : Return the total copy number value for the REGION passed as parameter in the region passed as parameter
    - getMinorCN : Return the minor copy number value for the REGION passed as parameter in the region passed as parameter
    - getLogR : Return the logR value for the REGION passed as parameter in the region passed as parameter
    - getPurity : Return the purity of the REGION passed as parameter
    - getPloidy : Return the ploidy of the REGION passed as parameter
    - getLikelyhood : Return the likelyhood/goodness of fit of the REGION passed as parameter
    - getGoodness : Return the likelyhood/goodness of fit of the REGION passed as parameter
"""

"""
REGION FORMAT:
{chr : [start, end, copy-number, total_copy_number, low_copy_number, logR_value]],
ploidy : float,
purity : float,
likelyhood : float}
"""

def getRegion(l) :
    """Return the start and end position from a chr list element in the REGION format

    Returns the first and second element from the list passed as parameter. That corresponds to start and end position in the chromsomal region passed as parameter

    Parameters:
        l (list) : Element in REGION["chr"] list

    Returns:
        list : The start and end position from the region in the format [start, end]
    """
    print l[0:1]
    return l[0:1]

def getCopyNumber(reg, chr, dc) :
    """Return the copy-number classification of the REGION passed as parameter, in the region and chromosome passed as parameter

    Returns the copy-number classification, that means (N)ormal, (A)mplification, (D)eletion, CNN-(L)OH, in the region passed as parameter. The function looks at the REGION passed as parameter

    Parameters:
        reg (list): region to check the copy number. Format of the list: [start, end]
        chr (string) : chromosome where the region is
        dc (REGION) : dict, in REGION format, to find the region, and return the copy number value

    Returns:
        char : 'N' when Normal copy number, 'A' when amplification, 'D' when deletion, 'L' when copy number neutral loss of Heterozygosity (CNN-LOH)
    """
    cnv = "N"
    if chr in dc.keys() :
        for r in dc[chr] :
            if r[0] <= reg[0] and r[1] >= reg[1] :
                cnv = r[2]
                break
    return cnv

def getTotalCN(reg, chr, dc) :
    """Return the total copy number value of the REGION passed as parameter, in the region and chromosome passed as parameter

    Returns the total copy number value in the region passed as parameter. The function looks at the REGION passed as parameter

    Parameters:
        reg (list): region to check the copy number. Format of the list: [start, end]
        chr (string) : chromosome where the region is
        dc (REGION) : dict, in REGION format, to find the region, and return the total copy number value

    Returns:
        int : Total copy number value in the region passed as parameter
    """
    cn = None
    if chr in dc.keys() :
        try :
            for r in dc[chr] :
                if r[0] <= reg[0] and r[1] >= reg[1] :
                    cn = r[3]
                    break
        except IndexError :
            cn = None
    return cn

def getMinorCN(reg, chr, dc) :
    """Return the minor copy number value of the REGION passed as parameter, in the region and chromosome passed as parameter

    Returns the minor copy number value in the region passed as parameter. The function looks at the REGION passed as parameter

    Parameters:
        reg (list): region to check the copy number. Format of the list: [start, end]
        chr (string) : chromosome where the region is
        dc (REGION) : dict, in REGION format, to find the region, and return the copy number value

    Returns:
        int : Minor copy number value in the region passed as parameter
    """
    cn = None
    if chr in dc.keys() :
        try :
            for r in dc[chr] :
                if r[0] <= reg[0] and r[1] >= reg[1] :
                    cn = r[4]
                    break
        except IndexError :
            cn = None
    return cn

def getLogR(reg, chr, dc) :
    """Return the logR value of the REGION passed as parameter, in the region and chromosome passed as parameter

    Returns the logR value in the region passed as parameter. The function looks at the REGION passed as parameter

    Parameters:
        reg (list): region to check the copy number. Format of the list: [start, end]
        chr (string) : chromosome where the region is
        dc (REGION) : dict, in REGION format, to find the region, and return the value

    Returns:
        float : LogR value in the region passed as parameter
    """
    logR = None
    if chr in dc.keys() :
        try :
            for r in dc[chr] :
                if r[0] <= reg[0] and r[1] >= reg[1] :
                    logR = r[5]
                    break
        except IndexError :
            logR = None
    return logR

def getPurity(dc) :
    """Get the purity of the REGION passed as parameter

    Returns:
        float : Purity value in the REGION passed as parameter
    """
    return dc["purity"]

def getPloidy(dc) :
    """Get the ploidy of the REGION passed as parameter

    Returns:
        float : Ploidy value of the region passed as parameter
    """
    return dc["ploidy"]

def getLikelyhood(dc) :
    """Get the likelyhood of the REGION passed as parameter

    Notice that this value is the same as the one obtained using getGoodness function

    Returns:
        float : log-likelyhood of the region passed as parameter
    """
    return dc["likelyhood"]

def getGoodness(dc) :
    """Get the goodnes of fit of the REGION passed as parameter

    Notice that this value is the same as the one obtained using getLikelyhood function

    Returns :
        float : Goodness of fit of the region passed as parameter
    """
    return dc["likelyhood"]

if __name__ == "__main__" :
    print "INFO:  This is a library file. Not able to be executed alone"
