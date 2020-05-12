#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MAIN: Programa principal. Muestra los pasos a seguir para poder hacer el analisis automatico de paneles en SMD. Usa la libreria creaLog para ejecutar el analisi
"""

import creaLog as cl
import subprocess

def GUI() :
    """
    Menu de interaccion con el usuario. Muestra las opciones de analisis disponibles
    """
    # TODO: Pintar-ho bonico
    # IDEA: Crear un menu amb els passos que es volen seguir per fer l'analisi mes personalitzat
    ruta = input("Introducir el path absoluto de la carpeta donde estan los FASTQ a analizar: ")
    cl.prepararPanel(ruta, ["copiar","fastqc", "aln", "recal", "bamqc", "coverage", "strelkaGerm", "vanno", "filtrar", "excel"])
    # TODO: Preguntar si executar la pipeline usant subprocess


if __name__ == "__main__" :
  GUI()
