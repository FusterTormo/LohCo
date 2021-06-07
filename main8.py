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
Check if each has vcf output (variant caller passed as parameter). If both samples (tumor and control) have vcf, check LOH reported by ASCAT2, FACETS, ascatNGS, PURPLE and Sequenza in the gene passed as parameter.
Do not do the same analysis in more than pair in the same submitter.
Report the statistics in a tsv with the columns
    * Submitter id
    * Analysis folder prefix (uuidTM_VS_uuidCN)
    * Germline worse variant found in the gene
    * Somatic worse variant found in the gene
    * Mean difference in VAF (Mean somaticVAF/germlineVAF in the variants reported by both vcf)
    * LOH categorization according to VAF mean difference
    * ASCAT2 LOH report in the gene
    * FACETS LOH report in the gene
    * ascatNGS LOH report in the gene
    * Sequenza LOH report in the gene
    * PURPLE LOH report in the gene
"""

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

# Functions
def getMaxMaf(ls) :
    """Return the maximum minor allele frequency from the list passed as parameter

    Calculates the maximum of the list passed as parameter. As "NA" or "." can be introduced, it checks if it is possible to convert each element to float

    Parameters
    ----------
        ls : list
            List of MAFs

    Returns
    -------
        float
            Maximum MAF. "NA" if it was not possible to get the maximum (e.g. All values in the list were "NA")
    """
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
    """Extract the variant allele frequency from the INFO column passed as parameter

    Parameters
    ----------
        info : str
            Column, from the vcf, where VAF (or reads forward/reverse) are stored
        vc : str
            Variant caller. This is used to know which schema must follow the info file

    Returns
    -------
        float
            VAF Calulated by dividing the number of alterated reads by the total of reads (in percentage)
    """
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
    """Calculate the mean vaf difference among the tumor and control variants passed as parameter

    Extracts the variants in common in tumor and in control and calculates the division in each variant. Later extracts the mean of all these variants

    Parameters
    ----------
        cn : dict
            Dict with all the variants found in control sample
        tm : dict
            Dict with all the variants found in tumor sample

    Returns
    -------
        float | str
            The mean of the cn/tm division in the variants common in tumor and control
            If there are not variants, returns "NA"
    """
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
    """Get the variants reported in the gene passed as parameter

    Using shell commands, gets the variants from the vcf file passed as parameter. Transforms the data in a list of dicts. The dict keys are: varType1 (exonic, intronic...),
    varType2 (nonsynonymous SNV, frameshift deletion...), maf (maximum MAF in all the MAFs reported), vaf (Variant Allele Frequency), GT (0/1, 1/1...)

    Parameters
    ----------
        path : str
            Path to vcf file. From this file the function will extract the variants in the gene passed as parameter. It will check if the file has Platypus format, or Strelka2 format
        gene : str
            Gene name. Variants from this gene will be extracted from the vcf passed as parameter

    Returns
    -------
        dict
            The dict has, chromosome-position as key; and variant type, maf, vaf and genotype as subkeys (see description for further details)
    """
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
    """Classify the sample according to the worst variant found

    Parameters
    ----------
        vars : dict
            List of variants reported in the patient
        maxMaf : float
            Maximum MAF to consider a nonsynonymous SNV variant as pathogenic. In case that maxMaf is passed as -1. All SNV variants will be considered as neutral

    Returns
    -------
        str
            + if the worst variant is considered pathogenic
            ? if the worst variant is a nonsynonymous SNV, and maximum MAF is higher than maxMaf
            - otherwise (the variants are considered not pathogenic)
    """
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
    """Find the ASCAT2 file in the ASCAT2 folder and get the LOH (or LOHs) in the region passed as parameter

    Check ASCAT2 LOH report in the region passed as parameter. If more than one ASCAT2 files are reported, it checks if all ASCAT2 files report the same aberration. In case the do not
    it reports "NF".

    Parameters
    ----------
        ascat : str
            ASCAT2 folder to get the ASCAT2 files
        reg : list
            List in REGION format to get check the region in which check the ASCAT2 output

    Returns
    -------
        str
            A, D, L, N If the same aberration is found in all ASCAT2 files: (A)mplification, (N)ormal, (L)OH, (D)eletion
            NF If no ASCAT2 files found, or if the aberration found in the region is not the same in all ASCAT2 files
    """
    cn = "NF"
    files = os.listdir(ascat)
    if len(files) == 1 :
        abs = "{}/{}".format(ascat, files[0])
        cn = getLOH(abs, "ascatarray", reg)
    else :
        for f in files :
            abs = "{}/{}".format(ascat, f)
            if cn == "" :
                cn = getLOH(abs, "ascatarray", reg)
            else :
                auxCn = getLOH(abs, "ascatarray", reg)
                if cn != auxCn : # If the ASCAT2 outputs does not output the same aberration, we do not include the result
                    cn =  "NF"
                    break
    return cn

def getLOH(path, program, gene) :
    """Find the copy number (or LOH) in the program passed as parameter in the gene passed as parameter too

    Additionally, get the purity in the sample calculated by the program

    Parameters
    ----------
        path : str
            Output file from the program. Here we will search the LOH (or CNV)
        program : str
            Program that has generated the output file passed in the previous parameter. Valid values are: facets, ascatngs, sequenza, purple, and ascatarray
        gene : list
            List that contains the chromosome, start and end position of the gene where we want to find the LOH

    Returns
    -------
        sol : str
            LOH found in the region. Values can be A, D, L, or N
        pur : float|str
            If the purity has been found in the output file, a float with the purity reported.
            Otherwise "NA"
    """
	sol = "Not found"
    pur = "Not found"
	if os.path.isfile(path):
		reg = lc.convert2region(path, program, "error")
        pur = lg.getPurity(reg)
		sol = lg.getCopyNumber(gene[1:3], gene[0], reg)

	return (sol, pur)

def doLoh(path, region) :
    """Get the program and the aberration reported by the program

    Given a file and a region, checks the LOH reported in the region passed as parameter. To do that, checks which program has been passed as parameter.

    Parameters
    ----------
        path : str
            File path to check the LOH in the region
        region : list
            List, in REGION, format to get the LOH reported

    Returns
    -------
        program : str
            LOH tool that has created the file passed as parameter. Possible values: facets, ascatngs, sequenza, purple, ascat2
        loh : str
            Aberration reported by the tool: (A)mplification, (L)OH, (D)eletion, (N)ormal, Not found (NF), Not Available (NA)
    """
    program = ""
    loh = ""
    purity = ""
    # Get from which program is the output
    if path.endswith("facets_comp_cncf.tsv") :
        program = "facets"
    elif path.endswith("copynumber.caveman.csv") :
        program = "ascatngs"
    elif path.endswith("_segments.txt") :
        program = "sequenza"
    elif path.endswith("TUMOR.purple.cnv.somatic.tsv") :
        program = "purple"
    else :
        program = "ascat2"
    if program == "ascat2" :
        loh = checkAscat(path, region)
        purity = "NA"
    else :
        loh, purity = getLOH(path, program, region)

    return (program, loh, purity)

def printRstring(var) :
    """Transform the dict to R schema

    Converts the dict, passed as parameter, that has the number of A, D, L, N, NF found to R syntax

    Parameters
    ----------
        var : dict
            Histogram that counts all the aberrations found

    Returns
    -------
        str
            Data in R format
    """
    str = ""
    for k in ["ascat2", "facets", "ascatngs", "sequenza"] :
        str += "{} <- c(".format(k)
        for v in ["A", "L", "D", "N", "NF"] :
            str += "{},".format(var[k][v])

        str = str.strip(",") + ")"
    return str

# Main program
def main(brcagene, genename, vcPath, maxMaf = 0.01) :
    """Main program

    Obtains, from the data base the cases (submitters) in OV cancer
    Searches in the database the pairs tumor-normal in each OV submitter
    For each pair tumor-normal, searches the number of vcf files available as well as the number of LOH output files available
    If there are 2 vcf files in the pair, gets the LOH reported by all the LOH tools. Additionally it calculates the difference in VAF between tumor and normal variants
    Prints all the data obtained in a tsv file
    Name of the file: {gene}_{maxMAF}_{variant_caller}_LOH.tsv
    """
    cont = 0
    data = []
    # Get submitters list
    with dbcon :
          cur = dbcon.cursor()
          q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
          cases = q.fetchall()

    print("INFO: Analysis will be done in {} cases".format(len(cases)))
    for c in cases :
        # Get the pairs: tumors and controls
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
                analysisdir = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0]) # The folder format for FACETS, ascatNGS, PURPLE and Sequenza is "[tumorUUID]_VS_[controlUUID]""
                vct = vcPath.format(tf)
                vcc = vcPath.format(cf)
                # Check if both variant calling files exist
                if os.path.isfile(vct) :
                    auxDc["vcfFiles"].append(vct)
                if os.path.isfile(vcc) :
                    auxDc["vcfFiles"].append(vcc)
                # Check the LOH outputs available
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
                purple = "{wd}/{folder}_PURPLE/TUMOR.purple.cnv.somatic.tsv".format(wd = workindir, folder = analysisdir)
                if os.path.isfile(purple) :
                    auxDc["lohFiles"].append(purple)
                # Get which pair tumor-control has more files available and save it for the next analyses
                if len(auxDc["vcfFiles"]) > len(submitter["vcfFiles"]) and len(auxDc["lohFiles"]) > len(submitter["vcfFiles"]) :
                    submitter = auxDc.copy()

        # Get the variants in the gene. If there are variant calling files for tumor and control sample
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
            # Classify the LOH depending on the mean obtained
            if meanVaf != "NA" :
                if meanVaf > 1.65 :
                    temp["vafVarCat"] = "L"
                else :
                    temp["vafVarCat"] = "N"
            else :
                temp["vafVarCat"] = "NF"
            # Classify the variants according to the pathogenicity
            temp["germVar"] = classifyVariants(cnVar, maxMaf)
            temp["somVar"] = classifyVariants(tmVar, maxMaf)
            # Get the LOH in all the available LOH files
            for fic in submitter["lohFiles"] :
                prog, loh, purity = doLoh(fic, brcagene)
                temp[prog] = loh
                temp["p_{}".format(prog)] = purity

            data.append(temp)
            print(temp)


    print("\nINFO: Final results. {} had a sample with both vcfs (tumor and control) done".format(cont))

    if vcPath.find("platypusGerm") > 0 :
        vc = "Platypus"
    elif vcPath.find("strelkaGerm") > 0 :
        vc = "Strelka2"

    filename = "{}_{}_{}_LOH.tsv".format(genename, maxMaf, vc)
    print("INFO: Writing output in {}".format(filename))
    with open(filename, "w") as fi :
        fi.write("submitter\tanalysis\tGermlineVar\tSomaticVar\tVAFvar\tLOHcat\tASCAT2\tFACETS\tascatNGS\tSequenza\tPURPLE\tp_FACETS\tp_AscatNGS\tp_Sequenza\tp_PURPLE\n")
        for l in data :
            fi.write("{sub}\t{anal}\t{germ}\t{som}\t{vaf}\t{catVaf}\t".format(sub = l["submitter"], anal = l["cmp"], germ = l["germVar"], som = l["somVar"], vaf = l["vafDif"], catVaf = l["vafVarCat"]))
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
                fi.write("{}\t".format(l["sequenza"]))
            else :
                fi.write("NA\t")
            if "purple" in l.keys() :
                fi.write("{}\t".format(l["purple"]))
            else :
                fi.write("NA\t")
            if "p_facets" in l.keys() :
                fi.write("{}\t".format(l["p_facets"]))
            else :
                fi.write("NA\t")
            if "p_ascatngs" in l.keys() :
                fi.write("{}\t".format(l["p_ascatngs"]))
            else :
                fi.write("NA\t")
            if "p_sequenza" in l.keys() :
                fi.write("{}\t".format(l["p_sequenza"]))
            else :
                fi.write("NA\t")
            if "p_purple" in l.keys() :
                fi.write("{}\n".format(l["p_purple"]))
            else :
                fi.write("NA\n")


if __name__ == "__main__" :
    brca1 = ["17", 43044295, 43170245]
    brca2 = ["13", 32315086, 32400266]
    maxMaf = 0.01
    variantCallingFile = "{}/platypusGerm/platypus.hg38_multianno.txt"
    #variantCallingFile = "{}/strelkaGerm/results/variants/strelka.hg38_multianno.txt"
    main(brca1, "BRCA1", variantCallingFile, maxMaf)
