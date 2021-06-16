#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Extracts the best Jaccard Index for each submitter"""

import sqlite3
import subprocess

def getJCC(sub, fil) :
    cmd = "grep {} {}".format(sub, fil)
    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = pr.communicate()
    out = std.decode().strip().split("\n")
    header = True
    col = 10 # Column number for Jaccard index
    dat = -1
    for l in out :
        aux = l.split("\t")
        try :
            dat = float(aux[col])
        except ValueError :
            pass
        except IndexError :
            print("ERROR: Not found data in {}. Grep found:".format(fil))
            print(l)

    if dat == -1 :
        dat = "NA"
    print("Returned: {} in {}".format(dat, fil))

    return dat

"""Main program"""

dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
# Get submitters list
with dbcon :
      cur = dbcon.cursor()
      q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
      cases = q.fetchall()

for c in cases :
    # For each submitter, extract the best JCC
    ascat = getJCC(c[0], "main6/ascat2VSascat2.tsv")
    dnacopy = getJCC(c[0], "main6/ascat2VSarray.tsv")
    facets = getJCC(c[0], "main6/ascat2VSfacets.tsv")
    ascatngs = getJCC(c[0], "main6/ascat2VSascatNGS.tsv")
    sequenza = getjcc(c[0], "main6/ascat2VSsequenza.tsv")
    purple = getJcc(c[0], "main6/ascat2VSpurple.tsv")
    break
