#!/usr/bin/env python

from Bio import SeqIO
import sys
import re

# silence Biopython warnings of improper indentation of Genbank files and CDS that are not multiples of 3
import warnings
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# define input file handle
genbankFile = open(sys.argv[1], 'r')

outputFile = open('output.gbk', 'w')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "ncRNA":
		
			if 'ydaO-yuaA' in ''.join(feature.qualifiers["product"]):
				feature.type = "misc_feature"
				noteString = ''.join(feature.qualifiers["product"]) + " element"
				feature.qualifiers["note"] = noteString
				feature.qualifiers.pop("product", None)
				locusTag = ''.join(feature.qualifiers["locus_tag"])
				print locusTag
				feature.qualifiers.pop("locus_tag", None)

	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()
