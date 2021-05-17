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

"""MAIN PROGRAM
Get the submitters from OV cancer, and the pairs tumor-control from each submitter.
Check if each has vcf output (variant caller passed as parameter). If both samples (tumor and control) have vcf, check LOH reported by ASCAT2, FACETS and Sequenza in the gene passed as parameter.
Do not do the same analysis in more than pair in the same submitter.
Report the statistics in a tsv with the columns
    * Submitter id
    * Analysis folder prefix (uuidTM_VS_uuidCN)
    * Germline worse variant found in the gene
    * Somatic worse variant found in the gene
    * Mean difference in VAF (Mean somaticVAF/germlineVAF in the variants reported by both vcf)
    * ASCAT2 LOH report in the gene
    * FACETS LOH report in the gene
    * ascatNGS LOH report in the gene
    * Sequenza LOH report in the gene
"""

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
            vaf = round(100 * num/den,2)
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
        mean = round(sum(common)/len(common), 2)
    else :
        mean = "NA"

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

def classifyVariants(vars, maxMaf) :
    classification = "-"
    for k, v in vars.items() :
        if v["varType1"] == "exonic" :
            if v["varType2"] in ctes.var_positive :
                classification = "+"
                break
            elif v["varType2"] in ctes.var_neutral : # Unknown variants will be reclassified according MAF passed as parameter (maxMaf)
                if v["maf"] == "NA" :
                    if maxMaf >= 0 :
                        classification = "+"
                        break
                    else : # When maxMaf is equal to -1, that means that mafs that are NA, will not be considered as pathogenic
                        classification = "?"
                else :
                    if v["maf"] <= maxMaf :
                        classification = "+"
                        break
                    else :
                        classification = "?"
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
    cont = 0
    data = []
    # Get submitters list
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
          cases = q.fetchall()

    print("INFO: Analysis will be done in {} cases".format(len(cases)))
    for c in cases :
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

        # Get the variants in the gene
        if len(submitter["vcfFiles"]) == 2 :
            cont += 1
            temp = {"submitter" : c[0], "cmp" : analysisdir}
            if cont % 100 == 0 :
                print("INFO: {} cases analysed".format(cont))
            tmVar = getVariant(submitter["vcfFiles"][0], genename)
            cnVar = getVariant(submitter["vcfFiles"][1], genename)
            # Get the VAF comparison to infer possible LOH
            meanVaf = getVafMean(cnVar, tmVar)
            temp["vafDif"] = meanVaf
            # Classify the variants according to the pathogenicity
            temp["germVar"] = classifyVariants(cnVar, maxMaf)
            temp["somVar"] = classifyVariants(tmVar, maxMaf)

            for fic in submitter["lohFiles"] :
                prog, loh = doLoh(fic, brcagene)
                temp[prog] = loh

            data.append(temp)


    print("\nINFO: Final results. {} had a sample with both vcfs (tumor and control) done".format(cont))
    print(data[-1])
    if vcPath.find("platypusGerm") > 0 :
        vc = "Platypus"
    elif vcPath.find("strelkaGerm") > 0 :
        vc = "Strelka2"
    filename = "{}_{}_{}_LOH.tsv".format(genename, maxMaf, vc)
    print("INFO: Writing output in {}".format(filename))
    with open(filename, "w") as fi :
        fi.write("submitter\tanalysis\tGermlineVar\tSomaticVar\tVAFvar\tASCAT2\tFACETS\tascatNGS\tSequenza\n")
        for l in data :
            fi.write("{sub}\t{anal}\t{germ}\t{som}\t{vaf}\t".format(sub = l["submitter"], anal = l["cmp"], germ = l["germVar"], som = l["somVar"], vaf = l["vafDif"]))
            if "ascat2" in l.keys() :
                fi.write("{}\t".format(l["ascat2"]))
            else :
                fi.write("NA\t")
            if "facets" in l.keys() :
                fi.write("{}\t".format(l["facets"]))
            else :
                fi.write("NA\t")
            if "ascatngs" in l.keys() :
                fi.write("{}\t".format(l["ascatngs"]))
            else :
                fi.write("NA\t")
            if "sequenza" in l.keys() :
                fi.write("{}\n".format(l["sequenza"]))
            else :
                fi.write("NA\n")
    # output = "INFO: Variables for R. Params: Gene: {}, MAF: {} ".format(genename, maxMaf)
    # if vcPath.find("platypusGerm") > 0 :
    #     output += "Variant caller: Platypus\n"
    # elif vcPath.find("strelkaGerm") > 0 :
    #     output += "Variant caller: Strelka2\n"
    # output += "\tPositive\n"
    # output += "\t\t{}\n".format(printRstring(dcPos))
    # output += "\tNegative\n"
    # output += "\t\t{}\n".format(printRstring(dcNeg))
    # output += "\tUnknown\n"
    # output += "\t\t{}\n".format(printRstring(dcNeu))
    # output += "\n------\n"
    # print(output)
    # sys.exit()

if __name__ == "__main__" :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    maxMaf = -1
    variantCallingFile = "{}/platypusGerm/platypus.hg38_multianno.txt"
    #variantCallingFile = "{}/strelkaGerm/results/variants/strelka.hg38_multianno.txt"
    main(brca1, "BRCA1", variantCallingFile, maxMaf)
