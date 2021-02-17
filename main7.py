#!/usr/bin/python
# -*- coding: utf-8 -*-

import multiprocessing as mlt
import os
import sqlite3
import subprocess
import sys

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

def getVariant(path, gene) :
    removableVars = ["intergenic", "intronic", "unkwnonw", "downstream", "upstream"]
    noMaf = []
    wMaf = {"varType1" : "", "varType2" : "", "maf" : 2}

    # Find the variants called in BRCA1
    cmd = "grep {gene} {vc}".format(vc = path, gene = gene)
    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = pr.communicate()
    out = std.decode().strip().split("\n")
    for tmp in out :
        aux = tmp.split("\t")
        try :
            varType = aux[5]
            varType2 = aux[8]
            if varType not in removableVars :
                maf = getMaxMaf(aux[10:39])
                if maf == "NA" :
                    noMaf.append({"varType1" : varType, "varType2" : varType2, "maf" : maf, "GT": aux[-1]})
                else :
                    if maf < wMaf["maf"] :
                        wMaf["varType1"] = varType
                        wMaf["varType2"] = varType2
                        wMaf["maf"] = maf
                        wMaf["GT"] = aux[-1]
        except IndexError :
            pass

    if wMaf["maf"] < 2 :
        return (wMaf, noMaf)
    else :
        return (None, noMaf)


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
def main(brcagene, genename, vcPath, maxMaf = 0.01, output) :
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
            gmmaf, gmnoMaf = getVariant(submitter["vcfFiles"][1], genename)
            # IDEA: Podria fer multiprocessing en aquesta busqueda
            #tmmaf, tmnoMaf = getVariant(submitter["vcfFiles"][0], "BRCA1")
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

    out = "INFO: Variables for R. Params: Gene: {}, MAF: {} ".format(genename, maxMaf)
    if vcPath.find("platypusGerm") > 0 :
        out += "Variant caller: Platypus\n"
    elif vcPath.find("strelkaGerm") > 0 :
        out += "Variant caller: Strelka2\n"
    out += "\tPositive\n"
    out += "\t\t{}\n".format(printRstring(dcPos))
    out += "\tNegative\n"
    out += "\t\t{}\n".format(printRstring(dcNeg))
    out += "\tUnknown\n"
    out += "\t\t{}\n".format(printRstring(dcNeu))


brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]
maxMaf = 0.05
variantCallingFile = "{}/platypusGerm/platypus.hg38_multianno.txt"
out1 = ""
out2 = ""
out3 = ""
#variantCallingFile = "{}/strelkaGerm/results/variants/strelka.hg38_multianno.txt"
p1 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.05, out1))
p2 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.03, out2))
p3 = mlt.Process(target=main, args = (brca1, "BRCA1", variantCallingFile, 0.0, out3))
p1.start()
p2.start()
p3.start()
#main(brca2, "BRCA2", variantCallingFile, maxMaf)
p1.join()
p2.join()
p3.join()
print(out1)
print(out2)
print(out3)
