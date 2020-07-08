#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of introduction of TCGA ASCAT2
"""

"""
The test includes:
    1. Open example outputs from sequenza, facets, ascatNGS, array, and ASCAT2
    2. Count the aberrations reported by each tool
    3. Get the copy number reported by the all outputs in BRCA1 and BRCA2 genes
"""

import libcomparison as lc

# Convert the files to REGION format
ascatngs = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs")
sequenza = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_Sequenza/TCGA-04-1332_segments.txt", "sequenza")
facets = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6c_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets")
array = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_E11_437726.grch38.seg.v2.txt", "array")
ascat = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/ASCAT2/TCGA-OV.79e63073-7d6d-456b-92c7-a3a7f0216ee7.ascat2.allelic_specific.seg.txt", "ascatarray")

print("INFO: Files opened successfully")
