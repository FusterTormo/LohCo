#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
MAIN: Funciones para modificar un nuevo manifest
"""

"""
AYUDA: Comandos para generar los manifest tal y como los necesitan los programas
    mv manifest.bed manifest_original.bed # Mantener una copia del original
    sort -k1,1 -k2,2n -V -s manifest_original.bed > manifestaux.bed # Ordenar el manifest por cromosoma y coordenada
    sed 's/chr//g' manifestaux.bed > manifest.bed # Eliminar los 'chr'
    bgzip -c manifest.bed > manifest.bed.gz # Comprimir el archivo para Strelka2. Puede que la ruta original de bgzip sea /opt/htslib-1.10.2/bin/bgzip
    tabix -p bed manifest.bed.gz # Crear un indice para el manifest. Puede que la ruta original de tabix sea /opt/htslib-1.10.2/bin/tabix
"""

def anotarManifest(ruta) :
    """Anotar un maniest. Anadiendo informacion del gen al que pertenece cada region

    Anotar un manifest pasado por parametro, usando la base de datos de UCSC (mysql) y el genoma de referencia hg19. Una vez anotado invoca a la funcion doListaGenes para crear el archivo
    gensAestudi.txt, necesario para calcular el coverage de cada gen en el panel

    Parameters
    ----------
        ruta : str
            Path donde esta el manifest que se quiere anotar
    """
    if os.path.isfile(ruta) :
        print("INFO: Anotando manifest")
        dir = os.path.dirname(ruta)
        cont = ""
        arx = "{}/manifest_anotado.bed".format(dir)
        with open(ruta, "r") as fi :
            for l in fi :
                aux = l.strip().split("\t")
                if len(aux) >= 3 :
                    if not aux[0].startswith("chr") :
                        aux[0] = "chr{}".format(aux[0])
                    # NOTA: En caso de querer cambiar el genoma de referencia, modificar en la siguiente linea hg19 por el genoma de interes
                    cmd = "mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -P 3306 -sN -D hg19 -e \"select name2 from refGene where chrom = '{chr}' AND txStart <= {nd} AND txEnd >= {st}\"".format(chr = aux[0], nd = aux[2], st = aux[1])
                    pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    std, err = pr.communicate()
                    out = std.decode()
                    genes = list(dict.fromkeys(out.strip().split("\n")))
                    cont += "{}\t{}\t{}\t{}\n".format(aux[0], aux[1], aux[2], ",".join(genes))
        with open(arx, "w") as fi :
            fi.write(cont)
            print("INFO: Manifest guardado como {}".format(arx))
            genes = "{}/gensAestudi.txt".format(dir)
            print("INFO: Creando el archivo gensAestudi.txt para calculos de coverage")
            doListaGenes(arx, genes)

def doListaGenes(manifest, genes) :
    """Crear el archivo con la lista de genes que hay dentro del manifest

    Este archivo se usa luego para calcular el coverage de las regiones de interes de cada gen dentro del manifest

    Parameters
    ----------
        manifest : str
            Ruta donde esta el manifest del cual se va a extraer la lista de genes
        genes : str
            Nombre y ruta que tendra el archivo con la lista de genes
    """
    listaGenes = []
    with open(manifest, "r") as fi :
        for l in fi :
            aux = l.split("\t")[3]
            aux = aux.strip()
            if aux not in listaGenes :
                listaGenes.append(aux)
    with open(genes, "w") as fi :
        fi.write("\n".join(listaGenes))
