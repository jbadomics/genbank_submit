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
	
		if feature.type == "CDS":
			
			# trimmedProduct = re.sub(" \(EC \d.*\)", "", str(feature.qualifiers["product"]).strip("['\"]"))
			# feature.qualifiers["product"] = trimmedProduct.replace('""', '"')
			# trimmedProduct = re.sub(" \(TC \d.*\)", "", str(feature.qualifiers["product"]).strip("['\"]"))
			# feature.qualifiers["product"] = trimmedProduct.replace('""', '"')
			product = ''.join(feature.qualifiers["product"]).replace(" [Streptomyces albus J1074]", "").replace(" [Streptomyces albus]", "").replace("[Streptomyces albus]", "")
			feature.qualifiers["product"] = product

			
	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()