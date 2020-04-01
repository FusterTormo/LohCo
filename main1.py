#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Llegir summary.md on hi ha les variants detectades en cadascuna de les mostres de OV i les mutacions que tenen en BRCA
Llegir els corresponents FACETS, ascatNGS i Sequenza
Comprovar que detecten en la regio. Si LOH o no
-----------------------------------------------
IDEA: Fer una taula amb tota la informacio.

	|Mostra|uuids[0:8]|Variant      |output FACETS        |output ascatNGS      |output Sequenza      |
	|-     |__VS__    |tumor|control|BRCA1|BRCA2|ATM|PALB2|BRCA1|BRCA2|ATM|PALB2|BRCA1|BRCA2|ATM|PALB2|
	|------|----------|-----|-------|-----|-----|---|-----|-----|-----|---|-----|-----|-----|---|-----|

? Fer un resum de la taula amb tota la informacio
    Positius FACETS_positius ascatNGS_positius Sequenza_positius
    Negatius FACETS_negatius ascatNGS_negatius Sequenza_negatius
"""

import sqlite3

# Vore l'estatus de les mostres i executar sequenza en aquelles que no s'haja executat
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
with dbcon :
	cur = dbcon.cursor()
	q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
	cases = q.fetchall()

print("INFO: {} cases found".format(len(cases)))
for c in cases :
	with dbcon :
		cur = dbcon.cursor()
		q = cur.execute("SELECT uuid FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
		tumors = q.fetchall()
		q = cur.execute("SELECT uuid FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
		controls = q.fetchall()
	print(c)
	print(tumors)
	print(controls)
	
