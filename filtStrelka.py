#!/usr/bin/python
# -*- coding: utf-8 -*-

import reAnnoFilt as anno

def calcularStrandBias(l) :
    """Calcular el numero total de reads en forward - numero total de reads en reverse

    Cuenta cuantos reads hay en forward y cuantos en reverse y los resta para que se pueda considerar si la variante tiene strand bias o no.

    Parameters
    ----------
        l : dict
            Variante en la que se quiere calcular el strand bias score
    """
    coef = "NP"
    auxF = l["ADF"].split(",")
    auxR = l["ADR"].split(",")
    forward = 0
    reverse = 0
    for a in auxF :
        forward += int(a)
    for a in auxR :
        reverse += int(a)
    coef = forward - reverse

    return str(coef)

def main(ruta, samplename = "noName") :
    """Programa principal. Filtra los datos de Strelka2 y guarda los datos en archivos .reanno.tsv

    Filtra los datos de Strelka2, previamente anotados con ANNOVAR. Anade informacion util a los datos ya existentes, como un score para strand bias, calidad del mapeo, VAF, coverage,
    enlaces para IGV, datos sobre la poblacion con MAF mas alta... Tambien llama a la funcion getWebInfo para recoger informacion de MAFs de myvariant.info. Finalmente, filtra los datos y
    guarda las variantes que no han pasado cada uno de los filtros en archivos .reannot.tsv. Estos archivos se usaran en data2excel para crear el excel final con todas las variantes y
    los filtros que han pasado. Tambien se crea un archivo de texto con estadisticas del numero de variantes que hay en cada proceso de filtrado

    Parameters
    ----------
        ruta : str
            Path donde esta el archivo hg19_multianno.txt (o hg38_multianno.txt) que se va a reanotar y filtrar.
        samplename : str, optional
            Nombre de la muestra. Este dato se incluira a los datos de interes
    """
    # Leer el archivo multianno de ANNOVAR y guardar los datos en un diccionario
    todas = anno.convertirData(ruta)
    # Agregar las columnas adicionales: "population_max", "predictor_summary". "Strand_bias_score", "Ref_depth", "Alt_depth", "VAF", "IGV_link", "sample"
    for p in todas :
        p["population_max"], p["population_max_name"] = anno.maximMaf(p)
        p["predictor_summary"] = anno.resumPredictors(p)
        p["Strand_bias_score"] = calcularStrandBias(p)
        p["Ref_depth"] = p["AD"].split(",")[0]
        p["Alt_depth"] = p["AD"].split(",")[1]
        if "DP" in p.keys() :
            p["VAF"] = float(p["Alt_depth"])/float((p["DP"])) * 100
        else :
            tmp = p["AD"].split(",")
            dp = 0
            for t in tmp :
                dp += int(t)
            p["VAF"] = float(p["Alt_depth"])/dp * 100
        p["IGV_link"] = "http://localhost:60151/goto?locus=chr{}:{}".format(p["Chr"], p["Start"])
        p["sample"] = samplename

    # Separar las variantes exonicas (consecuencia) y las splicing en un diccionario aparte
    conseq = anno.filtroConseq(todas)
    # Si el total de variantes es menor de 200 (numero arbitrario) re-anotar todas las variantes. En caso contrario, solo re-anotar las conseq
    if len(todas) <= 200 :
        print("INFO: Recogiendo informacion web para {} variantes".format(len(todas)))
        todas = anno.addWebInfo(todas)
    else :
        print("INFO: Recogiendo informacion web para {} variantes".format(len(conseq)))
        conseq = anno.addWebInfo(conseq)
    # Filtrar por MAF
    mafAlta, mafBaja = anno.filtrarMAF(conseq)
    # Filtrar por VAF
    vafAlta, vafBaja = anno.filtrarVAF(mafBaja)

    # Guardar los filtros en archivos de texto
    anno.guardarTabla(todas, "raw") # Todas las variantes recogidas en el vcf
    anno.guardarTabla(conseq, "conseq") # Variantes en regiones de splicing o exonicas (exluyendo sinonimas)
    anno.guardarTabla(mafAlta, "highMAF") # Variantes con una MAF>=0.01 en alguna columna de base de datos poblacional
    anno.guardarTabla(vafBaja, "lowVAF") # Variantes con una frecuencia alelica menor de 0.1
    anno.guardarTabla(vafAlta, "cand") # Variantes que han pasado todos los filtros anteriores

    # Mostrar por pantalla estadisticas basicas
    print("INFO: Variantes encontradas por Strelka2 {}".format(len(todas))) # Total de variantes reportadas por Strelka2
    print("INFO: Variantes conseq: {}".format(len(conseq)))
    print("INFO: MAF >= 0.1: {}".format(len(mafAlta)))
    print("INFO: MAF < 0.1: {}".format(len(mafBaja)))
    print("INFO: VAF baja: {}".format(len(vafBaja)))
    print("INFO: Variantes candidatas: {}".format(len(vafAlta)))

    # Agregar estadisticas de los tipos de variantes recogidos en el panel a la pestaña QC
    tiposVariantes = {}
    for t in todas :
        if t["Func.refGene"] == "exonic" :
            clave = "exonic_{}".format(t["ExonicFunc.refGene"])
        else :
            clave = t["Func.refGene"]
        if clave in tiposVariantes.keys() :
            tiposVariantes[clave] += 1
        else :
            tiposVariantes[clave] = 1

    with open("variants.stats.txt", "w") as fi :
        fi.write("{")
        fi.write("\'Totales\' : {},".format(len(todas)))
        fi.write("\'Conseq\' : {},".format(len(conseq)))
        fi.write("\'MAF_alta\' : {},".format(len(mafAlta)))
        fi.write("\'VAF_baja\' : {},".format(len(vafBaja)))
        fi.write("\'Candidatas\' : {}".format(len(vafAlta)))
        for k, v in tiposVariantes.items() :
            fi.write(",\'{}\' : {}".format(k, v))

        fi.write("}")

if __name__ == "__main__" :
    if len(sys.argv) > 2 :
        path = sys.argv[1]
        sample = sys.argv[2]
        main(path, sample)
    elif len(sys.argv) == 2 :
        main(path)
    else :
        print("ERROR: Numero de parametros invalido. Uso\n\tpython3 reAnnoFilt.py archivo_multianno.txt nombreMuestra")
