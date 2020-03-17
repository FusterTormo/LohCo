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
