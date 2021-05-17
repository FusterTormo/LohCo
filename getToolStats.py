#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Get, from a cancer and a tool defined by the user:
    * number of pairs tumor-control available
    * number of analyses succeed
    * number of analyses failed
"""

import sqlite3
import os

cancers = {"OV" : "/g/strcombio/fsupek_cancer2/TCGA_bam", "BRCA" : "/g/strcombio/fsupek_cancer3/TCGA_bam"}
tools = {
    "FACETS" : {"file" : "facets_comp_cncf.tsv", "suffix" : "_FACETS"},
    "ascatNGS" : {"file" : "copynumber.caveman.csv", "suffix" : "_ASCAT"},
    "Sequenza" : {"file" : "_segments.txt", "suffix" : "_Sequenza"},
    "PURPLE" : {"file" : "TUMOR.purple.cnv.gene.tsv", "suffix" : "_PURPLE"}
}

dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
file = input("Which tool you want to check the analyses?\nOptions: ({}) ".format(" - ".join(tools.keys())))
repo = input("Which cancer check the tool?\nOptions: ({}) ".format(" - ".join(cancers.keys())))

all = 0
done = 0
success = 0

with dbcon :
      cur = dbcon.cursor()
      q = cur.execute("SELECT submitter FROM patient WHERE cancer='{cancer}'".format(cancer = repo))
      cases = q.fetchall()

print("INFO: Total submitters in {} cancer: {}".format(repo, len(cases)))
if file in tools.keys() and repo in cancers.keys() :
    for c in cases :
        with dbcon :
            cur = dbcon.cursor()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
            tumors = q.fetchall()
            q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
            controls = q.fetchall()

        for tm in tumors :
            for cn in controls :
                all += 1
                prefix = "{}_VS_{}".format(tm[0].split("-")[0], cn[0].split("-")[0])
                wd = "{tmpath}/{cancer}/{submitter}/{prefix}{suffix}".format(tmpath = cancers[repo], cancer = repo, prefix = prefix, suffix = tools[file]["suffix"], submitter = c[0])
                if os.path.isdir(wd) :
                    done += 1
                    if file == "ascatNGS" or file == "Sequenza" :
                        aux = os.listdir(wd)
                        for a in aux :
                            if a.endswith(tools[file]["file"]) :
                                success += 1
                    elif file == "FACETS" or file == "PURPLE" :
                        aux = "{}/{}".format(wd, tools[file]["file"])
                        if os.path.isfile(aux) :
                            success += 1
    print("---------------\n--- Summary ---\n---------------\n")
    print("\t{} pairs tumor-control".format(all))
    print("\t{} analyses done".format(done))
    print("\t{} succeed, {} errors".format(success, done-success))


else :
    print("ERROR: Invalid option {} or {}".format(file, repo))
