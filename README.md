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
To get more examples, check the unitTest files:

* [unitTest1](../master/unitTest1.py)
* [unitTest2](../master/unitTest2.py)
* [unitTest3](../master/unitTest3.py)
* [unitTest4](../master/unitTest4.py)

Or to the main files:

* [main1](../master/main1.py)
