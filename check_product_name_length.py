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

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "CDS":
			
			locusTag = ''.join(feature.qualifiers["locus_tag"])
			product = ''.join(feature.qualifiers["product"])
			productStringLength = len(''.join(feature.qualifiers["product"]))
			if productStringLength > 100:
				print locusTag + "\t" + str(productStringLength) + "\t" + product

genbankFile.close()