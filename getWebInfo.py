import httplib2
import requests
import sys

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

def getMyVariant(chr, pos, ref, alt, na='NA') :
    mv = {'CADD_1000g_all' : na, 'CADD_1000g_afr' : na, 'CADD_1000g_amr' : na, 'CADD_1000g_eur' : na, 'CADD_1000g_eas' : na, 'CADD_1000g_sas' : na,
    'dbNSFP_1000g_all' : na, 'dbNSFP_1000g_afr' : na, 'dbNSFP_1000g_amr' : na, 'dbNSFP_1000g_eur' : na, 'dbNSFP_1000g_eas' : na, 'dbNSFP_1000g_sas' : na,
    'CADD_ESP6500_all' : na, 'CADD_ESP6500_ea' : na, 'CADD_ESP6500_aa' : na,
    'dbNSFP_esp6500_all' : na, 'dbNSFP_esp6500_ea' : na, 'dbNSFP_esp6500_aa' : na,
    'ExAC_ExAC_all' : na, 'ExAC_ExAC_afr' : na, 'ExAC_ExAC_amr' : na, 'ExAC_ExAC_eas' : na, 'ExAC_ExAC_fin' : na, 'ExAC_ExAC_nfe' : na, 'ExAC_ExAC_oth' : na, 'ExAC_ExAC_sas' : na,
    'dbNSFP_ExAC_all' : na, 'dbNSFP_ExAC_afr' : na, 'dbNSFP_ExAC_amr' : na, 'dbNSFP_ExAC_eas' : na, 'dbNSFP_ExAC_fin' : na, 'dbNSFP_ExAC_nfe' : na, 'dbNSFP_ExAC_oth' : na, 'dbNSFP_ExAC_sas' : na,
    'dbSNP_MAF' : na}

    if len(ref) == 1 and len(alt) == 1 and ref in ['A','C','G','T'] and alt in ['A','C','G','T'] :
        hgvs = chr + ':g.' + pos + ref + '>' + alt

        #Recoger datos de myvariant.info
        #Campos disponibles en myvariant.info http://docs.myvariant.info/en/latest/doc/data.html#available-fields
        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        campos = "cadd.1000g,dbnsfp.1000gp3,"
        campos += "cadd.esp,dbnsfp.esp6500,"
        campos += "exac,dbnsfp.exac,"
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
                                    mv['ExAC_all'] = "=%s/%s" % (a['exac']['ac']['ac'][it], a['exac']['an']['an'])
                                    break
                        else :
                            mv['ExAC_all'] = "=%s/%s" % (a['exac']['ac']['ac'], a['exac']['an']['an'])
                    if 'ac_afr' in a['exac']['ac'] and 'an_afr' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_afr'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_afr'] = "=%s/%s" % (a['exac']['ac']['ac_afr'][it], a['exac']['an']['an_afr'])
                                    break
                        else :
                            mv['ExAC_afr'] = "=%s/%s" % (a['exac']['ac']['ac_afr'], a['exac']['an']['an_afr'])
                    if 'ac_amr' in a['exac']['ac'] and 'an_amr' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_amr'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_amr'] = "=%s/%s" % (a['exac']['ac']['ac_amr'][it], a['exac']['an']['an_amr'])
                                    break
                        else :
                            mv['ExAC_amr'] = "=%s/%s" % (a['exac']['ac']['ac_amr'], a['exac']['an']['an_amr'])
                    if 'ac_eas' in a['exac']['ac'] and 'an_eas' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_eas'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_eas'] = "=%s/%s" % (a['exac']['ac']['ac_eas'][it], a['exac']['an']['an_eas'])
                                    break
                        else :
                            mv['ExAC_eas'] = "=%s/%s" % (a['exac']['ac']['ac_eas'], a['exac']['an']['an_eas'])
                    if 'ac_fin' in a['exac']['ac'] and 'an_fin' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_fin'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_fin'] = "=%s/%s" % (a['exac']['ac']['ac_fin'][it], a['exac']['an']['an_fin'])
                                    break
                        else :
                            mv['ExAC_fin'] = "=%s/%s" % (a['exac']['ac']['ac_fin'], a['exac']['an']['an_fin'])
                    if 'ac_nfe' in a['exac']['ac'] and 'an_nfe' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_nfe'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_nfe'] = "=%s/%s" % (a['exac']['ac']['ac_nfe'][it], a['exac']['an']['an_nfe'])
                                    break
                        else :
                            mv['ExAC_nfe'] = "=%s/%s" % (a['exac']['ac']['ac_nfe'], a['exac']['an']['an_nfe'])
                    if 'ac_sas' in a['exac']['ac'] and 'an_sas' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_sas'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_sas'] = "=%s/%s" % (a['exac']['ac']['ac_sas'][it], a['exac']['an']['an_sas'])
                                    break
                        else :
                            mv['ExAC_sas'] = "=%s/%s" % (a['exac']['ac']['ac_sas'], a['exac']['an']['an_sas'])
                    if 'ac_oth' in a['exac']['ac'] and 'an_oth' in a['exac']['an'] :
                        if isinstance(a['exac']['ac']['ac_oth'],(list, tuple)) :
                            for it in range(0,len(a['exac']['alleles'])) :
                                if (a['exac']['alleles'][it] == alt) :
                                    mv['ExAC_oth'] = "=%s/%s" % (a['exac']['ac']['ac_oth'][it], a['exac']['an']['an_oth'])
                                    break
                        else :
                            mv['ExAC_oth'] = "=%s/%s" % (a['exac']['ac']['ac_oth'], a['exac']['an']['an_oth'])

        else :
            print res
    return mv



def main() :
    header = False
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

if __name__ == "__main__" :
    main()
