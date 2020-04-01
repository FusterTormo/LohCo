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
import os

# Vore l'estatus de les mostres i executar sequenza en aquelles que no s'haja executat
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"
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
	for tm in tumors :
		for cn in controls :
			aux = tm[0].split("-")[0]
			aux2 = cn[0].split("-")[0]
			folder = "{wd}/{sub}/{tumor}_VS_{control}".format(wd = wd, sub = c[0], tumor = aux, control = aux2)
			# Comprovar l'status de l'analisi de Sequenza
			aux = "{}_Sequenza".format(folder)
			print("Comprovant {}".format(aux))
			if os.path.isfolder("{}_Sequenza".format(folder)) :
				# Si existeix la carpeta es que s'ha fet el seqz de la mostra. Comprovar si s'ha fet el sequenza
				aux = "{}_Sequenza".format(folder)
				if not os.path.isfile("{}/{}_segments.txt".format(aux, c[0])) :
					aux = "{}/{}_segments.txt".format(aux, c[0])
					print("Executar plotSequenza en {}".format(aux))
					break
			else :
				print("{} no existeix".format(aux))
