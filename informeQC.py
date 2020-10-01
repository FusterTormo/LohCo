#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import shutil

pathTemplate = "/home/ffuster/AUP/informes/informe.html"
pathCss = "/home/ffuster/AUP/informes/estils.css"

# Comprobar que las carpetas existen
if os.path.isdir("fastqc") and os.path.isdir("coverage") :
    dirsFQ = [] # Lista con los directorios donde se han guardado los datos de FastQC
    covs = [] # Lista con los graficos de coverage de cada gen
    mostra = os.getcwd().split("/")[-1] # Identificador de la muestra
    # Recoger el nombre de los directorios del output de FastQC
    for root, dirsFQ, files in os.walk("fastqc") :
        break

    # Recoger los png que del coverage
    for root, dirs, files in os.walk("coverage") :
        for f in files :
            if f.endswith("png") and f != "coverage.png" :
                covs.append(f)
        break
    # Ordenar los png alfabeticamente
    covs.sort()
    dirsFQ.sort()

    # Crear el contenido que tendra la seccion de estadisticas generales
    tab1 = "<table>"
    tab2 = "<table>"
    tab3 = "<table>"
    if os.path.isfile("alnQC.txt") :
        tab1 += "<thead><tr><th>Reads en</th><th>Read no.</th></thead>\n<tbody>\n"
        with open("alnQC.txt", "r") as fi :
            aux = fi.read()
        aln = eval(aux)
        tab1 += "<tr><th>FASTQ</th><td>{}</td></tr>\n".format(aln["FASTQ"])
        tab1 += "<tr><th>BAM</th><td>{} ({:.2f} %)</td></tr>\n".format(aln["BAM"], 100*float(aln["BAM"])/float(aln["FASTQ"]))
        tab1 += "<tr><th>ON target</th><td>{} ({:.2f} %)</td></tr>\n".format(aln["ON"], 100*float(aln["ON"])/float(aln["BAM"]))
        tab1 += "<tr><th>OFF target</th><td>{} ({:.2f} %)</td></tr>\n".format(aln["OFF"], 100*float(aln["OFF"])/float(aln["BAM"]))
        tab1 += "</tbody>\n"
    tab1 += "</table></div>\n"
    if os.path.isfile("coverage.txt") :
        tab2 += "<thead><tr><th colspan='2' style='text-align: center'>Coverage general</th></thead>\n<tbody>\n"
        with open("coverage.txt", "r") as fi :
            aux = fi.read()
        cv = eval(aux)
        tab2 += "<tr><th>M&iacute;nim</th><td>{}</td></tr>".format(cv["minimo"])
        tab2 += "<tr><th>M&aacute;xim</th><td>{}</td></tr>".format(cv["maximo"])
        tab2 += "<tr><th>Mediana</th><td>{}</td></tr>".format(cv["mediana"])
        tab2 += "<tr><th>Mitjana</th><td>{:.2f}</td></tr>".format(float(cv["media"]))
        tab2 += "<tr><th>% Bases cov = 0</th><td>{:.2f}</td></tr>".format(cv["bases0"])
        tab2 += "<tr><th>% Bases cov &le; 30</th><td>{:.2f}</td></tr>".format(cv["bases30"])
        tab2 += "<tr><th>% Bases cov &le; 100</th><td>{:.2f}</td></tr>".format(cv["bases100"])
        tab2 += "<tr><th>% Bases cov &le; 500</th><td>{:.2f}</td></tr>".format(cv["bases500"])
        tab2 += "<tr><th>% Bases cov &le; 1000</th><td>{:.2f}</td></tr>".format(cv["bases1000"])
        tab2 += "</tbody>"
    tab2 += "</table>\n"
    if os.path.isfile("coverageGeneStats.tsv") :
        header = False
        with open("coverageGeneStats.tsv", "r") as fi:
            for l in fi :
                aux = l.strip().split("\t")
                if not header :
                    tab3 += "<thead><tr><th>Gen</th><th>M&iacute;nim</th><th>M&agrave;xim</th><th>Mitjana</th><th>Mediana</th></tr></thead><tbody>\n"
                    header = True
                else :
                    tab3 += "<tr><td>{gen}</td><td>{min}</td><td>{max}</td><td>{mean:.2f}</td><td>{md}</td></tr>\n".format(
                        gen = aux[0].strip("\""), min = aux[1], max = aux[2], mean = float(aux[3]), md = aux[4])
            tab3 += "</tbody>"
    tab3 += "</table></div>\n"

    # Crear el contenido de la seccion que contiene los coverages de los genes
    covsTxt = ""
    for c in covs :
        covsTxt += "\t\t<a href=\"coverage/{png}\" title=\"Clica a l'imatge per ampliar-la\">\n\t\t\t<img src=\"coverage/{png}\">\n\t\t</a>\n".format(png = c)

    with open(pathTemplate, "r") as fi :
        txt = fi.read()

    # Guardar les dades del informe en un HTML
    with open("informe.html", "w") as fi :
        fi.write(txt.format(idMostra = mostra, dirFASTQ1 = dirsFQ[0], dirFASTQ2 = dirsFQ[1], imgCoverageGens = covsTxt, taulesStats = tab1+tab2+tab3))

    # Copiar el full d'estils en la carpeta
    shutil.copyfile(pathCss, "estils.css")
    print("INFO: Creado informe de calidad de coverage y FASTQ")

else :
    print("ERROR: No encontradas las carpetas necesarias para hacer el informe: {dir}/fastqc y {dir}/coverage".format(dir = os.getcwd()))
