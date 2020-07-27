# Comparing different tools for LOH analysis

## Examples
### Test 1) Comparison between FACETS and CNV array obtained from TCGA

This tests includes

1. Read output example from array and FACETS
2. Fragment the regions to get a common group of regions
3. Calculate 4x4 table
4. Extract confusion matrices for (A)mplification, (D)eletion, and (N)ormal copy number

| Array/FACETS  | Amplification | LOH | Normal | Deletion |
| ------------- | ------------- | --- | ------ | -------- |
| Amplification | ------------- | --- | ------ | -------- |
| LOH           | ------------- | --- | ------ | -------- |
| Normal        | ------------- | --- | ------ | -------- |
| Deletion      | ------------- | --- | ------ | -------- |

```python
import libcomparison as compi
import libstatistics as sts

print("TEST 1) Extract data")
ar = compi.convert2region("input_examples/arrays/TCGA-13-0887_tumor.txt", "array")
fa = compi.convert2region("input_examples/facets_comp_cncf.tsv", "FACETS")
print("TEST 2) Search regions in common for the study")
regs = compi.getFragments(ar, fa)
print("TEST 3) Create the 4x4 comparison table")
dc = compi.doComparison(regs, ar, fa)
print(sts.printTable(dc, "Array", "FACETS", False))
print("TEST 4) Get the confusion matrix of each aberration. TCGA arrays cannot detect LOH")
c1, c2 = sts.calculateCounts(dc)
dicContingency = sts.doContingency(dc, ["A", "D", "N"])
print("\tAmplificacion\n\t{}\n\n".format(dicContingency["A"]))
print("\tDelecion\n\t{}\n\n".format(dicContingency["D"]))
print("\tNormal\n\t{}\n".format(dicContingency["N"]))
print("TEST 5) Jaccard index de cada una de las aberraciones")
jci2 = compi.doComparison2(regs, ar, fa)
dicJaccard = sts.jaccardIndex(jci2, ["A", "D", "N"])
print("\tAmplificacion\n\t{}\n\n".format(dicJaccard["A"]))
print("\tDelecion\n\t{}\n\n".format(dicJaccard["D"]))
print("\tNormal\n\t{}\n".format(dicJaccard["N"]))
```

### Test 2) Comparison between FACETS and ascatNGS output

This test includes

1. Read FACETS and ascatNGS example outputs
2. Fragment the regions to get a common group of regions
3. Get ploidy and purity from each output
4. Plot the correlation between logR
5. ggplot of total copy number
6. ggplot of minor copy number
7. Write a bed file with the regions got in each file and the common regions

```python
import libcomparison as compi
import libstatistics as sts
import libgetters as gt
import os

print("TEST 1) Extract data")
ascat = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/TCGA-09-0369/TCGA-09-0369_40e311a4_VS_f4441d6e/H_GP-09-0369-01A-01W-0372-09-1.copynumber.caveman.csv", "ascat")
facets = compi.convert2region("/home/labs/solelab/ffuster2/Desktop/doctorat/cas_estudi/input_examples/TCGA-09-0369/TCGA-09-0369_40e311a4_VS_f4441d6e/facets_comp_cncf.tsv", "FACETS")
print("TEST 2) Divide the regions to get common regions")
regs = compi.getFragments(facets, ascat)
print("TEST 3) Get the ploidy")
print("\tFACETS: {}".format(gt.getPloidy(facets))
print("\tascatNGS: {}".format(gt.getPloidy(ascat)))
print("TEST 3) Get the purity")
print("\tFACETS: {}".format(gt.getPurity(facets))
print("\tascatNGS: {}".format(gt.getPurity(ascat)))
print("TEST 4) Plot logR concordance")
try :
    sts.logRcomp(regs, facets, ascat, "FACETS", "ASCAT")
except ValueError :
    print("ERROR: Cannot create the logR plot")
print("TEST 5) Plot copy number counts using ggplot")
sts.doGGplotFiles(facets, ascat, "FACETS", "ASCAT")
print("TEST 6) Write a bed file with all the reported regions and the common regions")
sts.print2Bed(ascat, "ascatNGS", facets, "FACETS", regs)
```
### Test 5) Comparison between ASCAT2, SNP-array, FACETS, ascatNGS, and Sequenza in a particular region

This test opens the data from this 5 different sources, counts the number of regions reported by each aberration in each program. Finally obtains the copy number reported by each tool in BRCA1 and BRCA2 genes, respectively


```python
import libcomparison as lc
import libgetters as lg
import libstatistics as ls

# BRCA1/2 gene coordinates
brca1 = ["17", 43044295, 43170245]
brca2 = ["13", 32315086, 32400266]

# Convert the files to REGION format
ascatngs = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs", "error")
sequenza = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_Sequenza/TCGA-04-1332_segments.txt", "sequenza", "error")
facets = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets", "error")
array = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_E11_437726.grch38.seg.v2.txt", "array")
ascat = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/ASCAT2/TCGA-OV.79e63073-7d6d-456b-92c7-a3a7f0216ee7.ascat2.allelic_specific.seg.txt", "ascatarray")

print("INFO: Files opened successfully")
print("INFO: Number of aberrations reported by each tool")
print("Array: {}".format(ls.countsXtool(array)))
print("ASCAT2: {}".format(ls.countsXtool(ascat)))
print("ascatNGS: {}".format(ls.countsXtool(ascatngs)))
print("FACETS: {}".format(ls.countsXtool(facets)))
print("Sequenza: {}".format(ls.countsXtool(sequenza)))
print("----------------------------------------")
print("INFO: Copy number reported in BRCA1")
print("ASCAT2  : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascat)))
print("Array   : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], array)))
print("ascatNGS: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], ascatngs)))
print("FACETS  : {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], facets)))
print("Sequenza: {}".format(lg.getCopyNumber(brca1[1:3], brca1[0], sequenza)))
print("----------------------------------------")
print("INFO: Copy number reported in BRCA2")
print("ASCAT2  : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascat)))
print("Array   : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], array)))
print("ascatNGS: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], ascatngs)))
print("FACETS  : {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], facets)))
print("Sequenza: {}".format(lg.getCopyNumber(brca2[1:3], brca2[0], sequenza)))

```


To get more examples, check the unitTest files:

* [unitTest1](../master/unitTest1.py)
* [unitTest2](../master/unitTest2.py)
* [unitTest3](../master/unitTest3.py)
* [unitTest4](../master/unitTest4.py)
* [unitTest5](../master/unitTest5.py)

Or to the main files:

* [main1](../master/main1.py)
* [main2](../master/main2.py)
* [main3](../master/main3.py)
* [main4](../master/main4.py)
* [main5](../master/main5.py)
