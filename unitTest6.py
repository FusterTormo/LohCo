#!/usr/bin/python
# -*- coding: utf-8 -*-

import libcomparison as lc
import libgetters as lg
import libstatistics as ls

# Convert the files to REGION format
# ascatngs = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_ASCAT/H_GP-04-1332-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs", "error")
# sequenza = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_Sequenza/TCGA-04-1332_segments.txt", "sequenza", "error")
facets = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/90cf56c6_VS_f4b549d0_FACETS/facets_comp_cncf.tsv", "facets", "error")
# array = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/Array/QUANT_p_TCGA_Batch12_AFFX_GenomeWideSNP_6_E11_437726.grch38.seg.v2.txt", "array")
# ascat = lc.convert2region("/g/strcombio/fsupek_cancer2/TCGA_bam/OV/TCGA-04-1332/ASCAT2/TCGA-OV.79e63073-7d6d-456b-92c7-a3a7f0216ee7.ascat2.allelic_specific.seg.txt", "ascatarray")

# Count the FACETS aberrations
print(ls.countsXtool(facets))
fg = lc.getFragments(facets, facets)
print(fg)
