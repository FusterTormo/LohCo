#!/usr/bin/python
# -*- coding: utf-8 -*-

import multiprocessing as mlt
import os
import sqlite3
import subprocess
import sys
import time

import main1 as lib
import libconstants as ctes

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

# Functions
def getMaxMaf(ls) :
    maf = -1
    for l in ls :
        try :
            aux = float(l)
            if aux >= maf :
                maf = aux
        except ValueError :
            pass
    if maf == -1 :
        maf = "NA"

    return maf

def getVaf(info, vc) :
    num = -1
    den = -1
    vaf = -1
    if vc == "platypus" :
        for tmp in info.split(";") :
            key, val = tmp.split("=")
            if key == "TR" :
                try :
                    num = float(val)
                except ValueError : # To fix multivariant sites. Get one of the multivariant
                    aux = val.split(",")[0]
                    num = float(aux)
            if key == "TC" :
                den = float(val)
        if num >= 0 and den > 0 :
            vaf = 100 * num/den
    else :
        print("WARNING: Strelka not implemented yet")
        sys.exit()

    return vaf

def getVafMean(cn, tm) :
    common = []
    aux = 0
    for c in cn.keys() :
        if c in tm.keys() :
            aux = tm[c]["vaf"]/cn[c]["vaf"]
            common.append(aux)
    if len(common) > 0 :
        mean = sum(common)/len(common)
    else :
        mean = "NA"
    print("{} variants in common. {} mean vaf".format(len(common), mean))
    return mean

def getVariant(path, gene) :
    removableVars = ["intergenic", "intronic", "unkwnonw", "downstream", "upstream"]
    noMaf = []
    var = {}

    # Find the variants called in the gene passed as parameter
    cmd = "grep {gene} {vc}".format(vc = path, gene = gene)
    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = pr.communicate()
    out = std.decode().strip().split("\n")
    for tmp in out :
        aux = tmp.split("\t")
        try :
            chrom = aux[0]
            pos = aux[1]
            key = "{}-{}".format(chrom, pos)
            varType = aux[5]
            varType2 = aux[8]
            if varType not in removableVars :
                maf = getMaxMaf(aux[10:39])
                if path.endswith("platypus.hg38_multianno.txt") :
                    vc = "platypus"
                elif path.endswith("strelka.hg38_multianno.txt") :
                    vc = "strelka2"

                vaf = getVaf(aux[47], vc)
                var[key] = {"varType1" : varType, "varType2" : varType2, "maf" : maf, "vaf" : vaf, "GT" : aux[-1]}
        except IndexError :
            pass

    return var

def classifyVariants(maf, noMaf, maxMaf) :
    classification = "-"
    if maf is not None :
        if maf["varType1"] == "exonic" :
            if maf["varType2"] in ctes.var_negative :
                classification = "-"
            elif maf["varType2"] in ctes.var_positive :
                classification = "+"
            else :
                if maf["maf"] <= maxMaf :
                    classification = "+"
                else :
                    classification = "?"
        elif maf["varType1"] == "splicing" :
            classification = "+"

    if maxMaf < 0 : # If the MAF has a value below 0, we consider the variants with no MAF as unknown
        if len(noMaf) > 0 :
            classification = "?"
    else :
        if classification == "-" or classification == "?" :
            for v in noMaf :
                if v["varType1"] == "exonic" :
                    if v["varType2"] in ctes.var_neutral or v["varType2"] in ctes.var_positive :
                        classification = "+"
                        break
                elif v["varType1"] == "splicing" :
                    classification = "+"
                    break
    return classification

def checkAscat(ascat, reg) :
    cn = "NF"
    files = os.listdir(ascat)
    if len(files) == 1 :
        abs = "{}/{}".format(ascat, files[0])
        cn = lib.getLOH(abs, "ascatarray", reg)
    else :
        for f in files :
            abs = "{}/{}".format(ascat, f)
            if cn == "" :
                cn = lib.getLOH(abs, "ascatarray", reg)
            else :
                auxCn = lib.getLOH(abs, "ascatarray", reg)
                if cn != auxCn : # If the ASCAT2 outpus does not output the same aberration, we do not include the result
                    cn =  "NF"
                    break
    return cn

def doLoh(path, region) :
    program = ""
    loh = ""
    # Get from which program is the output
    if path.endswith("facets_comp_cncf.tsv") :
        program = "facets"
    elif path.endswith("copynumber.caveman.csv") :
        program = "ascatngs"
    elif path.endswith("_segments.txt") :
        program = "sequenza"
    else :
        program = "ascat2"
    if program == "ascat2" :
        loh = checkAscat(path, region)
    else :
        loh = lib.getLOH(path, program, region)

    return (program, loh)

def printRstring(var) :
    str = ""
    for k in ["ascat2", "facets", "ascatngs", "sequenza"] :
        str += "{} <- c(".format(k)
        for v in ["A", "L", "D", "N", "NF"] :
            str += "{},".format(var[k][v])

        str = str.strip(",") + ")"
    return str

# Main program
def main(brcagene, genename, vcPath, maxMaf = 0.01) :
    totalPos = 0
    dcPos = {"ascat2" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "facets" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0},
    "ascatngs" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "sequenza" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}}
    totalNeg = 0
    dcNeg = {"ascat2" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "facets" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0},
    "ascatngs" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "sequenza" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}}
    totalNeu = 0
    dcNeu = {"ascat2" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "facets" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0},
    "ascatngs" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}, "sequenza" : {"L" : 0, "A" : 0, "D" : 0, "N" : 0, "NF" : 0}}
    # Get submitters list
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
          cases = q.fetchall()

    print("INFO: Analysis done in {} cases".format(len(cases)))
    for c in cases :
        if (totalPos + totalNeg + totalNeu) % 100 == 0 :
            print("{} cases analysed".format(totalPos + totalNeg + totalNeu))

        with dbcon :
            cur = dbcon.cursor()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
            tumors = q.fetchall()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
            controls = q.fetchall()

        # Get the pair tumor/control with more files
        submitter = {"tumor" : "", "control" : "", "vcfFiles" : [], "lohFiles" : []}
        for tm in tumors :
            for cn in controls :
                auxDc = {"tumor" : tm[0], "control" : cn[0], "vcfFiles" : [], "lohFiles" : []}
                tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
                cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])
                workindir = "{wd}/{sub}".format(wd = wd, sub = c[0])
                analysisdir = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, and Sequenza is "[tumorUUID]_VS_[controlUUID]""
                vct = vcPath.format(tf)
                vcc = vcPath.format(cf)
                # Check both variant calling files exist
                if os.path.isfile(vct) :
                    auxDc["vcfFiles"].append(vct)
                if os.path.isfile(vcc) :
                    auxDc["vcfFiles"].append(vcc)
                folder = "{}/ASCAT2".format(workindir)
                if os.path.isdir(folder) and len(os.listdir(folder)) > 0:
                    auxDc["lohFiles"].append(folder)
                facets = "{wd}/{folder}_FACETS/facets_comp_cncf.tsv".format(wd = workindir, folder = analysisdir)
                if os.path.isfile(facets) :
                    auxDc["lohFiles"].append(facets)
                ascatngs = lib.findAscatName("{wd}/{folder}_ASCAT/".format(wd = workindir, folder = analysisdir))
                if ascatngs != "Not found" :
                    auxDc["lohFiles"].append(ascatngs)
                sequenza = "{wd}/{folder}_Sequenza/{case}_segments.txt".format(folder = analysisdir, case = c[0], wd = workindir)
                if os.path.isfile(sequenza) :
                    auxDc["lohFiles"].append(sequenza)
                if len(auxDc["vcfFiles"]) > len(submitter["vcfFiles"]) and len(auxDc["lohFiles"]) > len(submitter["vcfFiles"]) :
                    submitter = auxDc.copy()

        # Get the variants in BRCA1 and BRCA2 genes
        if len(submitter["vcfFiles"]) == 2 :
            cnVar = getVariant(submitter["vcfFiles"][1], genename)
            tmVar = getVariant(submitter["vcfFiles"][0], genename)
            getVafMean(cnVar, tmVar)

            continue
            # IDEA: Podria fer multiprocessing en aquesta busqueda
            # Classify the variants according to the pathogenicity
            varClass = classifyVariants(gmmaf, gmnoMaf, maxMaf)
            # Get LOH in the region
            try :
                prog1, loh1 = doLoh(submitter["lohFiles"][0], brcagene)
                if varClass == "+" :
                    dcPos[prog1][loh1] += 1
                elif varClass == "-" :
                    dcNeg[prog1][loh1] += 1
                elif varClass == "?" :
                    dcNeu[prog1][loh1] += 1
            except IndexError :
                prog1 = None
                loh1 = None
            try :
                prog2, loh2 = doLoh(submitter["lohFiles"][1], brcagene)
                if varClass == "+" :
                    dcPos[prog2][loh2] += 1
                elif varClass == "-" :
                    dcNeg[prog2][loh2] += 1
                elif varClass == "?" :
                    dcNeu[prog2][loh2] += 1
            except IndexError :
                prog2 = None
                loh2 = None
            try :
                prog3, loh3 = doLoh(submitter["lohFiles"][2], brcagene)
                if varClass == "+" :
                    dcPos[prog3][loh3] += 1
                elif varClass == "-" :
                    dcNeg[prog3][loh3] += 1
                elif varClass == "?" :
                    dcNeu[prog3][loh3] += 1
            except IndexError :
                prog3 = None
                loh3 = None
            try :
                prog4, loh4 = doLoh(submitter["lohFiles"][3], brcagene)
                if varClass == "+" :
                    dcPos[prog4][loh4] += 1
                elif varClass == "-" :
                    dcNeg[prog4][loh4] += 1
                elif varClass == "?" :
                    dcNeu[prog4][loh4] += 1
            except IndexError :
                prog4 = None
                loh4 = None
            # Count the number of each
            if varClass == "+" :
                totalPos += 1
            elif varClass == "-" :
                totalNeg += 1
            elif varClass == "?" :
                totalNeu += 1
            # # # IDEA: Podria fer multiprocessing en la busqueda de cada arxiu LOH
            # print("{} -> {} - {}".format(submitter["lohFiles"][0], prog1, loh1))
            # print("{} -> {} - {}".format(submitter["lohFiles"][1], prog2, loh2))
            # print("{} -> {} - {}".format(submitter["lohFiles"][2], prog3, loh3))
            # print(dcPos)

    print("\nINFO: Final results. {} were able to analyse".format(totalNeg + totalNeu + totalPos))
    if vcPath.find("platypusGerm") > 0 :
        print("INFO: Using Platypus as variant caller")
    elif vcPath.find("strelkaGerm") > 0 :
        print("INFO: Using Strelka2 as variant caller")
    print("{} cases considered positive in {}".format(totalPos, genename))
    for k in dcPos.keys() :
        print("\t{} ({} analyses)".format(k, sum(dcPos[k].values())))
        nf = totalPos - sum(dcPos[k].values()) + dcPos[k]["NF"]
        dcPos[k]["NF"] = nf

        for key, value in dcPos[k].items() :
            print("\t\t{} -> {} found".format(key, value))
        print("\t\t\t{:.2f}% LOH".format(100 * (dcPos[k]["D"] + dcPos[k]["L"])/totalPos))

    print("\n{} cases considered negative in {}".format(totalNeg, genename))
    for k in dcNeg.keys() :
        print("\t{} ({} analyses)".format(k, sum(dcNeg[k].values())))
        nf = totalNeg - sum(dcNeg[k].values()) + dcNeg[k]["NF"]
        dcNeg[k]["NF"] = nf

        for key, value in dcNeg[k].items() :
            print("\t\t{} -> {} found".format(key, value))
        print("\t\t\t{:.2f}% LOH".format(100 * (dcNeg[k]["D"] + dcNeg[k]["L"])/totalNeg))

    print("\n{} cases considered unknown in {}".format(totalNeu, genename))
    for k in dcNeu.keys() :
        print("\t{} ({} analyses)".format(k, sum(dcNeu[k].values())))
        nf = totalNeu - sum(dcNeu[k].values()) + dcNeu[k]["NF"]
        dcNeu[k]["NF"] = nf

        for key, value in dcNeu[k].items() :
            print("\t\t{} -> {} found".format(key, value))
        print("\t\t\t{:.2f}% LOH".format(100 * (dcNeu[k]["D"] + dcNeu[k]["L"])/totalNeu))

    output = "INFO: Variables for R. Params: Gene: {}, MAF: {} ".format(genename, maxMaf)
    if vcPath.find("platypusGerm") > 0 :
        output += "Variant caller: Platypus\n"
    elif vcPath.find("strelkaGerm") > 0 :
        output += "Variant caller: Strelka2\n"
    output += "\tPositive\n"
    output += "\t\t{}\n".format(printRstring(dcPos))
    output += "\tNegative\n"
    output += "\t\t{}\n".format(printRstring(dcNeg))
    output += "\tUnknown\n"
    output += "\t\t{}\n".format(printRstring(dcNeu))
    output += "\n------\n"
    print(output)

if __name__ == "__main__" :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    maxMaf = 0.05
    variantCallingFile = "{}/platypusGerm/platypus.hg38_multianno.txt"
    #variantCallingFile = "{}/strelkaGerm/results/variants/strelka.hg38_multianno.txt"
    # p1 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.05))
    # p2 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.03))
    # p3 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.0))
    # p4 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, -1))
    # p1.start()
    # time.sleep(40)
    # p2.start()
    # time.sleep(40)
    # p3.start()
    # time.sleep(40)
    # p4.start()
    main(brca1, "BRCA1", variantCallingFile, maxMaf)
