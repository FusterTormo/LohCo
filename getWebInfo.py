#!/usr/bin/python
# -*- coding: utf-8 -*-

import httplib2
import requests
import sys

"""
WARNING - Funcion desechada temporalmente. La informacion del VEP por ahora no es necesaria y puede que este script este obsoleto
def getVep() :
    #Recoger datos de ENSEMBL VEP (solo datos de splicing)
    #http://grch37.rest.ensembl.org/vep/human/hgvs/ZRSR2:c.827+1G>A?content-type=application/json&MaxEntScan=true&GeneSplicer=true&dbscSNV=true
    #La variante corresponde al paciente CD34-775 secuenciado en la tanda11

    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/human/hgvs/ZRSR2:c.827+1G>A?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, params={"MaxEntScan" : True, "GeneSplicer" : True, "dbscSNV" : True})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    print repr(decoded)
"""

def getMyVariant(chr, pos, ref, alt, na='NA') :
    """
    Recoger los datos desde myvariant.info para una mutacion puntual dada

    Consulta en myvariant.info la informacion disponible de base de datos poblacionales de la variante que se ha pasado por parametro. Recoge los datos de
    1000genomes (desde CADD y dbNSFP)
    EVS (desde CADD y dbNSFP)
    ExAC (desde ExAC y dbNSFP)
    gNOMAD (genome y exome)
    dbSNP

    Parameters
    ----------
        chr : str
            Nombre del cromosoma. Se valida/corrige en caso de que el formato no sea correcto para myvariant.info
        pos : str
            Posicion donde esta la variante que se quiere consultar
        ref : str
            Base de referencia que ha sido modificada. Se comprueba que sea un unico caracter (A, C, T o G)
        alt : str
            Base observada en lugar de la de referencia. Se comprueba que sea un unico caracter (A, C, T o G)
        na : str, optional
            Texto que aparecera cuando no se encuentra informacion en un campo. Por defecto 'NA'
    Returns
    -------
        dict
            Diccionario con todos los datos recogidos para la variable. Los nombres de los indices se pueden ver en la variable siguiente
    """
    mv = {'CADD_1000g_all' : na, 'CADD_1000g_afr' : na, 'CADD_1000g_amr' : na, 'CADD_1000g_eur' : na, 'CADD_1000g_eas' : na, 'CADD_1000g_sas' : na,
    'dbNSFP_1000g_all' : na, 'dbNSFP_1000g_afr' : na, 'dbNSFP_1000g_amr' : na, 'dbNSFP_1000g_eur' : na, 'dbNSFP_1000g_eas' : na, 'dbNSFP_1000g_sas' : na,
    'CADD_ESP6500_all' : na, 'CADD_ESP6500_ea' : na, 'CADD_ESP6500_aa' : na,'dbNSFP_esp6500_all' : na, 'dbNSFP_esp6500_ea' : na, 'dbNSFP_esp6500_aa' : na,
    'ExAC_ExAC_all' : na, 'ExAC_ExAC_afr' : na,'ExAC_ExAC_amr' : na, 'ExAC_ExAC_eas' : na, 'ExAC_ExAC_fin' : na, 'ExAC_ExAC_nfe' : na, 'ExAC_ExAC_oth' : na, 'ExAC_ExAC_sas' : na,
    'dbNSFP_ExAC_all' : na, 'dbNSFP_ExAC_afr' : na, 'dbNSFP_ExAC_amr' : na, 'dbNSFP_ExAC_eas' : na, 'dbNSFP_ExAC_fin' : na, 'dbNSFP_ExAC_nfe' : na, 'dbNSFP_ExAC_oth' : na,
    'dbNSFP_ExAC_sas' : na,
    'gNOMAD_Exome_all' : na, 'gNOMAD_Exome_afr' : na, 'gNOMAD_Exome_amr' : na, 'gNOMAD_Exome_asj' : na, 'gNOMAD_Exome_eas' : na, 'gNOMAD_Exome_fin' : na, 'gNOMAD_Exome_nfe' : na,
    'gNOMAD_Exome_oth' : na, 'gNOMAD_Exome_popmax' : na, 'gNOMAD_Exome_raw' : na, 'gNOMAD_Exome_sas' : na,
    'gNOMAD_Genome_all' : na, 'gNOMAD_Genome_afr'  : na, 'gNOMAD_Genome_amr' : na, 'gNOMAD_Genome_asj' : na, 'gNOMAD_Genome_eas' : na, 'gNOMAD_Genome_fin' : na, 'gNOMAD_Genome_nfe' : na,
    'gNOMAD_Genome_oth' : na, 'gNOMAD_Genome_popmax' : na, 'gNOMAD_Genome_raw' : na,
    'dbSNP_MAF' : na}

    # Comrpobar si el cromosoma tiene el formato adecuado
    if not chr.startswith('chr') :
        chr = "chr{}".format(chr)

    if len(ref) == 1 and len(alt) == 1 and ref in ['A','C','G','T'] and alt in ['A','C','G','T'] :
        hgvs = chr + ':g.' + pos + ref + '>' + alt

        #Recoger datos de myvariant.info
        #Campos disponibles en myvariant.info http://docs.myvariant.info/en/latest/doc/data.html#available-fields
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        campos = "cadd.1000g,dbnsfp.1000gp3,"
        campos += "cadd.esp,dbnsfp.esp6500,"
        campos += "exac,dbnsfp.exac,"
        campos += "gnomad_exome,gnomad_genome,"
        campos += "dbsnp.gmaf"
        #cadd.esp,dbnsfp.esp6500,dbnsfp.1000gp3,dbnsfp.exac,dbsnp.rsid,dbsnp.gmaf,exac,exac_nontcga,snpeff"
        params = 'ids=' + hgvs + '&fields=%s' % (campos)
        res, con = h.request('http://myvariant.info/v1/variant', 'POST', params, headers=headers)
        #Comprobar respuesta correcta

        if res.status == 200 :
            try :
                a = eval(con)[0] #Convertir la respuesta en un array. Cada elemento del array sera un diccionario
            except NameError :
                a = []
            if 'dbsnp' in a :
                if 'gmaf' in a['dbsnp'] :
                    mv['dbSNP_MAF'] = str(a['dbsnp']['gmaf'])
            if 'cadd' in a :
                if '1000g' in a['cadd'] :
                    if 'af' in a['cadd']['1000g'] :
                        mv['CADD_1000g_all'] = str(a['cadd']['1000g']['af'])
                    if 'afr' in a['cadd']['1000g'] :
                        mv['CADD_1000g_afr'] = str(a['cadd']['1000g']['afr'])
                    if 'amr' in a['cadd']['1000g'] :
                        mv['CADD_1000g_amr'] = str(a['cadd']['1000g']['amr'])
                    if 'eur' in a['cadd']['1000g'] :
                        mv['CADD_1000g_eur'] = str(a['cadd']['1000g']['eur'])
                    if 'eas' in a['cadd']['1000g'] :
                        mv['CADD_1000g_eas'] = str(a['cadd']['1000g']['eas'])
                    if 'sas' in a['cadd']['1000g'] :
                        mv['CADD_1000g_sas'] = str(a['cadd']['1000g']['sas'])
                if 'esp' in a :
                    if 'af' in a['cadd']['esp']:
                        mv['CADD_ESP6500_all'] = str(a['cadd']['esp']['af'])
                    if 'eur' in a['cadd']['esp'] :
                        mv['CADD_ESP6500_ea'] = str(a['cadd']['esp']['eur'])
                    if 'afr' in a['cadd']['esp'] :
                        mv['CADD_ESP6500_aa'] = str(a['cadd']['esp']['afr'])
            if 'dbnsfp' in a :
                if '1000gp3' in a['dbnsfp'] :
                    if 'af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_all'] = str(a['dbnsfp']['1000gp3']['af'])
                    if 'afr_af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_afr'] = str(a['dbnsfp']['1000gp3']['afr_af'])
                    if 'amr_af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_amr'] = str(a['dbnsfp']['1000gp3']['amr_af'])
                    if 'eur_af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_eur'] = str(a['dbnsfp']['1000gp3']['eur_af'])
                    if 'eas_af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_eas'] = str(a['dbnsfp']['1000gp3']['eas_af'])
                    if 'sas_af' in a['dbnsfp']['1000gp3'] :
                        mv['dbNSFP_1000g_sas'] = str(a['dbnsfp']['1000gp3']['sas_af'])
                if 'esp6500' in a['dbnsfp'] :
                    if 'ea_af' in a['dbnsfp']['esp6500'] :
                        mv['dbNSFP_esp6500_ea'] = str(a['dbnsfp']['esp6500']['ea_af'])
                    if 'aa_af' in a['dbnsfp']['esp6500'] :
                        mv['dbNSFP_esp6500_aa'] = str(a['dbnsfp']['esp6500']['aa_af'])
                    #Por ahora no hay resultados de esp6500_all en dbNSFP
                if 'exac' in a['dbnsfp'] :
                    if 'af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_all'] = str(a['dbnsfp']['exac']['af'])
                    if 'afr_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_afr'] = str(a['dbnsfp']['exac']['afr_af'])
                    if 'amr_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_amr'] = str(a['dbnsfp']['exac']['amr_af'])
                    if 'eas_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_eas'] = str(a['dbnsfp']['exac']['eas_af'])
                    if 'fin_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_fin'] = str(a['dbnsfp']['exac']['fin_af'])
                    if 'nfe_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_nfe'] = str(a['dbnsfp']['exac']['nfe_af'])
                    #No hay datos para la poblacion Other en dbNSFP actualmente
                    if 'sas_af' in a['dbnsfp']['exac'] :
                        mv['dbNSFP_ExAC_sas'] = str(a['dbnsfp']['exac']['sas_af'])
            if 'gnomad_exome' in a :
                if 'af' in a['gnomad_exome'] :
                    if 'af' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_all'] = a['gnomad_exome']['af']['af']
                    if 'af_afr' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_afr'] = a['gnomad_exome']['af']['af_afr']
                    if 'af_amr' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_amr'] = a['gnomad_exome']['af']['af_amr']
                    if 'af_asj' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_asj'] = a['gnomad_exome']['af']['af_asj']
                    if 'af_eas' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_eas'] = a['gnomad_exome']['af']['af_eas']
                    if 'af_fin' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_fin'] = a['gnomad_exome']['af']['af_fin']
                    if 'af_nfe' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_nfe'] = a['gnomad_exome']['af']['af_nfe']
                    if 'af_oth' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_oth'] = a['gnomad_exome']['af']['af_oth']
                    if 'af_popmax' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_popmax'] = a['gnomad_exome']['af']['af_popmax']
                    if 'af_raw' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_raw'] = a['gnomad_exome']['af']['af_raw']
                    if 'af_sas' in a['gnomad_exome']['af'] :
                        mv['gNOMAD_Exome_sas'] = a['gnomad_exome']['af']['af_sas']
            if 'gnomad_genome' in a :
                if 'af' in a['gnomad_genome'] :
                    if 'af' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_all'] = a['gnomad_genome']['af']['af']
                    if 'af_afr' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_afr'] = a['gnomad_genome']['af']['af_afr']
                    if 'af_amr' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_amr'] = a['gnomad_genome']['af']['af_amr']
                    if 'af_asj' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_asj'] = a['gnomad_genome']['af']['af_asj']
                    if 'af_eas' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_eas'] = a['gnomad_genome']['af']['af_eas']
                    if 'af_fin' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_fin'] = a['gnomad_genome']['af']['af_fin']
                    if 'af_nfe' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_nfe'] = a['gnomad_genome']['af']['af_nfe']
                    if 'af_oth' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_oth'] = a['gnomad_genome']['af']['af_oth']
                    if 'af_popmax' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_popmax'] = a['gnomad_genome']['af']['af_popmax']
                    if 'af_raw' in a['gnomad_genome']['af'] :
                        mv['gNOMAD_Genome_raw'] = a['gnomad_genome']['af']['af_raw']
            if 'exac' in a :
                #La base de datos de exac no devuelve frecuencias alelicas, sino numero de veces que el alelo aparece en sus bases de datos y numero de veces del alelo normal
                #Si en la posicion hay variantes multialelicas, exac devuelve el contaje de todos los alelos. Hay que comprobar cual es el alelo que se busca
                #Para sacar la frecuencia hay que dividir el campo ac (allele count) entre el campo an (allele number)
                if 'ac' in a['exac'] and 'an' in a['exac']:
                    if 'ac' in a['exac']['ac'] and 'an' in a['exac']['an'] :
                        #Guarda los datos en formato formula de Excel, para asi no perder los datos de contaje alelico
                        if isinstance(a['exac']['ac']['ac'],(list, tuple)) : #ExAC ha devuelto una lista de alelos, buscamos el que corresponde a nuestra variante y lo guardamos
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac'][it])/float(a['exac']['an']['an'])
                                    mv['ExAC_ExAC_all'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac'])/float(a['exac']['an']['an'])
                            mv['ExAC_ExAC_all'] = "{:.2f}".format(aux)
                    if 'ac_afr' in a['exac']['ac'] and 'an_afr' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_afr'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_afr'][it])/float(a['exac']['an']['an_afr'])
                                    mv['ExAC_ExAC_afr'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_afr'])/float(a['exac']['an']['an_afr'])
                            mv['ExAC_ExAC_afr'] = "{:.2f}".format(aux)
                    if 'ac_amr' in a['exac']['ac'] and 'an_amr' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_amr'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_amr'][it])/float(a['exac']['an']['an_amr'])
                                    mv['ExAC_ExAC_amr'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_amr'])/float(a['exac']['an']['an_amr'])
                            mv['ExAC_ExAC_amr'] = "{:.2f}".format(aux)
                    if 'ac_eas' in a['exac']['ac'] and 'an_eas' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_eas'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_eas'][it])/float(a['exac']['an']['an_eas'])
                                    mv['ExAC_ExAC_eas'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_eas'])/float(a['exac']['an']['an_eas'])
                            mv['ExAC_ExAC_eas'] = "{:.2f}".format(aux)
                    if 'ac_fin' in a['exac']['ac'] and 'an_fin' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_fin'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_fin'][it])/float()
                                    mv['ExAC_ExAC_fin'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_fin'])/float(a['exac']['an']['an_fin'])
                            mv['ExAC_ExAC_fin'] = "{:.2f}".format(aux)
                    if 'ac_nfe' in a['exac']['ac'] and 'an_nfe' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_nfe'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_nfe'][it])/float(a['exac']['an']['an_nfe'])
                                    mv['ExAC_ExAC_nfe'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_nfe'])/float(a['exac']['an']['an_nfe'])
                            mv['ExAC_ExAC_nfe'] = "{:.2f}".format(aux)
                    if 'ac_sas' in a['exac']['ac'] and 'an_sas' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_sas'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_sas'][it])/float(a['exac']['an']['an_sas'])
                                    mv['ExAC_ExAC_sas'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_sas'])/float(a['exac']['an']['an_sas'])
                            mv['ExAC_ExAC_sas'] = "{:.2f}".format(aux)
                    if 'ac_oth' in a['exac']['ac'] and 'an_oth' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_oth'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    aux = 100*float(a['exac']['ac']['ac_oth'][it])/float(a['exac']['an']['an_oth'])
                                    mv['ExAC_ExAC_oth'] = "{:.2f}".format(aux)
                                    break
                        else :
                            aux = 100*float(a['exac']['ac']['ac_oth'])/float(a['exac']['an']['an_oth'])
                            mv['ExAC_ExAC_oth'] = "{:.2f}".format(aux)

        else :
            print(res)
    return mv



def main() :
    # Dummy unitary test
    header = False
    print(getMyVariant("6", "26093141", "26093141", "G", "A"))
    """
    with open("MO739/filtro4.allInfo","r") as fi :
        for l in fi :
            if not header :
                header = True
                continue
            else :
                e = l.split("\t")
                if e[8] == "nonsynonymous SNV" :
                    crom = e[0]
                    start = e[1]
                    ref = e[3]
                    alt = e[4]
                    val = getMyVariant(crom, start, ref, alt)
    """

if __name__ == "__main__" :
    main()
