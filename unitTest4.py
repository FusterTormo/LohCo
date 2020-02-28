#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Example of introduction of sequenza
"""

import libcomparison as lc

# sequenza = lc.convert2region("input_examples/sequenzaExample/TCGA-04-1343_segments.txt", "sequenza")
# facets = lc.convert2region("input_examples/TCGA-04-1343/TCGA-04-1343_9ece2e42_VS_65188600/facets_comp_cncf.tsv", "facets")
# ascat = lc.convert2region("input_examples/TCGA-04-1343/TCGA-04-1343_9ece2e42_VS_65188600/H_GP-04-1343-01A-01W-0488-09-1.copynumber.caveman.csv", "ascatngs")
print("INFO: Arxius oberts satisfactoriament")

lc.getFragments(None, None, [["2", 12341, 1234123474], ["2", 43234, 44444]])
