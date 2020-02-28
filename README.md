# Comparing different tools for LOH analysis

## Examples
### Test 1)

1. Read output example from array and FACETS.
2. Fragment the regions to get a common group of regions.
3. Calculate 4x4 table
|Array/FACETS |Amplification|LOH|Normal|Deletion|
|-------------|-------------|---|------|--------|
|Amplification|-------------|---|------|--------|
|LOH          |-------------|---|------|--------|
|Normal       |-------------|---|------|--------|
|Deletion     |-------------|---|------|--------|
4. Extract confusion matrices for (A)mplification, (D)eletion, and (N)ormal copy number

```
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
