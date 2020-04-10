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
import sys
import subprocess

# Constants
dbcon = sqlite3.connect("/g/strcombio/fsupek_cancer2/TCGA_bam/info/info.db")
wd = "/g/strcombio/fsupek_cancer2/TCGA_bam/OV"

# Vore l'estatus de les mostres i executar sequenza en aquelles que no s'haja executat
def checkSequenza() :
	sequenza = "/home/ffuster/Scripts/plotSequenza.R"
	cont = 0
	errors = 0
	messages = []
	with dbcon :
		cur = dbcon.cursor()
		q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
		cases = q.fetchall()

	print("INFO: {} cases found".format(len(cases)))
	for c in cases :
		cont += 1
		# Recollir la informacio dels bams i el sexe que te el cas registrats
		with dbcon :
			cur = dbcon.cursor()
			q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
			tumors = q.fetchall()
			q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
			controls = q.fetchall()
			q = cur.execute("SELECT gender FROM patient WHERE submitter='{}'".format(c[0]))
			gender = q.fetchone()[0]
		for tm in tumors :
			for cn in controls :
				aux = tm[0].split("-")[0]
				aux2 = cn[0].split("-")[0]
				folder = "{wd}/{sub}/{tumor}_VS_{control}".format(wd = wd, sub = c[0], tumor = aux, control = aux2)
				seq = "{}_Sequenza".format(folder)
				# Comprovar l'status de l'analisi de Sequenza
				if os.path.isdir(seq) :
					# Si existeix la carpeta es que s'ha fet el seqz de la mostra. Comprovar si s'ha fet el sequenza
					aux = "{}_Sequenza".format(folder)
					if not os.path.exists("{}/{}_segments.txt".format(aux, c[0])) :
						cmd = "Rscript {seqScript} {cas} {folder} {gender}".format(cas = c[0], folder = aux, gender = gender, seqScript = sequenza)
						print("INFO: Checked {cases} cases. Running {com}".format(cases = cont, com = cmd))
						proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
						out, err = proc.communicate()
						if proc.returncode != 0 :
							errors += 1
							messages.append(cmd)
							print("ERROR: While runing {}\nDescription:\n{}".format(cmd, err))
				else :
					aux = "{wd}/{sub}/{uuid}/{bam}".format(wd = wd, sub = c[0], uuid = tm[0], bam = tm[1])
					aux2 = "{wd}/{sub}/{uuid}/{bam}".format(wd = wd, sub = c[0], uuid = cn[0], bam = cn[1])
					cmd = "/home/ffuster/Scripts/runSequenza.sh {tumor} {control} {dir}".format(dir = seq, tumor = aux, control = aux2)
					print("INFO: Executing {}".format(cmd))
					proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
					out, err = proc.communicate()
					if proc.returncode != 0 :
						errors += 1
						messages.append(cmd)
						print("ERROR: While running {}\nDescription:\n{}".format(cmd, err))
					else :
						cmd = "Rscript {seqScript} {cas} {folder} {gender}".format(cas = c[0], folder = seq, gender = gender, seqScript = sequenza)
						proc = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
						out, err = proc.communicate()
						if proc.returncode != 0 :
							errors += 1
							messages.append(cmd)
							print("ERROR: While runing {}\nDescription:\n{}".format(cmd, err))

	print("INFO: Found {} errors in total.\nCommands failed:".format(errors))
	for m in messages :
		print("-> {}".format(m))

	print("Total errors: {}".format(errors))

# Buscar la pitjor variant reportada en el gen passat per parametre
def getWorst(vcf, gene) :
	found = False
	classifier = ["synonymous SNV", "nonsynonymous SNV", "nonframeshift substitution", "nonframeshift deletion",  "frameshift deletion", "frameshift insertion", "stopgain"]
	level = -1
	order = []
	worst = ""
	with open(vcf, "r") as fi :
		for l in fi :
			aux = l.split("\t") # Get the gene name
			print(l)
			if aux[6] == gene :
				found = True
				if classifier.index(aux[8]) > level :
					level = classifier.index(aux[8])
					worst = l
			elif found :
				break
	return worst

# Preparar la taula amb les mostres, les seues variants en els gens a estudi (BRCA1, BRCA2, ATM i PALB2) i el Copy number associat en cadascuna de les regions
def prepareTable() :
	cont = 0
	with dbcon :
		cur = dbcon.cursor()
		q = cur.execute("SELECT submitter FROM patient WHERE cancer='OV'")
		cases = q.fetchall()
	for c in cases :
		cont += 1
		print("INFO: {} cases checked".format(cont))
		# Recollir la informacio dels bams i el sexe que te el cas registrats
		with dbcon :
			cur = dbcon.cursor()
			q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Tumor%'".format(c[0]))
			tumors = q.fetchall()
			q = cur.execute("SELECT uuid, bamName FROM sample WHERE submitter='{}' AND tumor LIKE '%Normal%'".format(c[0]))
			controls = q.fetchall()
		for tm in tumors :
			for cn in controls :
				print("Checking Case {}. Tumor id {}. Control id {}".format(c[0], tm[0], cn[0]))
				tf = "{wd}/{sub}/{tumor}".format(wd = wd, sub = c[0], tumor = tm[0])
				platypus = "{}/platypusGerm/platypus.hg38_multianno.txt".format(tf)
				variant = getWorst(platypus, "ATM") # TODO: Fer per tots els gens
				cf = "{wd}/{sub}/{control}".format(wd = wd, sub = c[0], control = cn[0])

prepareTable()