#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
MAIN: Extrae parametros de calidad de los archivos vcf. Estos parametros son:
* Ka/Ks ratio: Mide el ratio entre mutaciones neutrales (sinonimas) y selectivas (no sinonimas) (https://en.wikipedia.org/wiki/Ka/Ks_ratio)
* Reads en forward que contienen una variante vs Reads en reverse con la variante. De acuerdo a https://pubmed.ncbi.nlm.nih.gov/28209900/ sirve para comprobar si ha habido errores en la preparacion del panel
* Ratios de cambios de base en una tabla para poder crear un histograma. Los cambios de bases se simplifican a 6: C>A, C>G, C>T, T>A, T>C, y T>G
"""

import re
import sys

def getRatios(vcf, imprimir = True) :
    """Calcular las ratios nombradas en el comentario anterior usando el archivo .hg19_multianno.txt crudo creado con ANNOVAR

    Parameters
    ----------
        vcf : str
            Ruta donde esta el archivo vcf en el que quiere evaluar la calidad
        imprimir : bool
            Mostrar por pantalla los resultados de los ratios obtenidos

    Returns
    -------
        float
            Ratio Ka/Ks
        int
            Diferencia entre los reads forward y los reverse que contienen variantes
        dict
            Diccionario con el numero de veces que se ha encontrado cada cambio de base
    """
    cmp = {"A" : "T", "C" : "G", "T" : "A", "G" : "C"} # Diccionario con las bases complementarias
    iref = None
    ialt = None
    isnv = None
    ihead = None
    ival = None
    isb = None
    hist = {"C>A" : 0, "C>G" : 0, "C>T" : 0, "T>A" : 0, "T>G" : 0, "T>C" : 0} # Histograma con el numero de cambios de bases
    # Numerador y denominador para el ratio dNS/dS
    dns = 0
    ds = 0
    # Numerador y denominador para el ratio FW/RV
    fw = 0
    rv = 0
    with open(vcf, "r") as fi :
        for l in fi :
            aux = l.strip().split("\t")
            if iref == None : # Obtener los numeros de columna que tienen los datos que se necesitan
                iref = aux.index("Ref")
                ialt = aux.index("Alt")
                isnv = aux.index("ExonicFunc.refGene")
                ihead = aux.index("Otherinfo9")
                ival = aux.index("Otherinfo10")
            else :
                if aux[isnv] == "nonsynonymous SNV" or aux[isnv] == "synonymous SNV" :
                    # Contar el numero de variantes sinonimas y no sinonimas hay en el vcf
                    if aux[isnv] == "nonsynonymous SNV" :
                        dns += 1
                    elif aux[isnv] == "synonymous SNV" :
                        ds += 1

                    # Contar el numero de cambios
                    if aux[iref] == "A" or aux[iref] == "G" : # En estas bases cuenta su complementario por lo mencionado en Alexandrov et al. (2013)
                        head = "{}>{}".format(cmp[aux[iref]], cmp[aux[ialt]])
                    else :
                        head = "{}>{}".format(aux[iref], aux[ialt])
                    hist[head] += 1

                    # Contar los reads con variantes forward y los reads con variantes en reverse
                    aux2 = aux[ihead].split(":") # Separar la cabecera de la columna FORMAT del archivo vcf
                    isb = aux2.index("SB")
                    aux2 = aux[ival].split(":")
                    sb = aux2[isb].split(",") # Orden de los valores en la columna 'SB' de GATK: referenceFW, referenceRV, alteratedFW, alteratedRV
                    fw += int(sb[2])
                    rv += int(sb[3])
    try :
        dnds = float(dns)/float(ds)
    except ZeroDivisionError :
        dnds = "NA"
        
    fwrv = int(fw) - int(rv)

    if imprimir :
        print("\nParametros de calidad de {}".format(vcf))
        print("------------------------------------------------")
        print("Ratio variantes no sinonimas/sinonimas:\t{:.2f}".format(dnds))
        print("Ratio variantes en Forward y Reverse:\t{}".format(fwrv))
        print("Histograma con los cambios SNV: ")
        for k,v in hist.items() :
            print("\t{}\t{}".format(k, v))

    return (dnds, fwrv, hist)

if __name__ == "__main__" :
    if len(sys.argv) > 1 :
        getRatios(sys.argv[1])
    else :
        print("\nINFO: Extraer parametros de calidad de un vcf creado con GATK Mutect2 y anotado con ANNOVAR\n\nUSO:\n\t./vcfQC.py archivo.hg19_multianno.txt\n\n")
