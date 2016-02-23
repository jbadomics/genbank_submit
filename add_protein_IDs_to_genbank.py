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
			
			if "note" in feature.qualifiers:
				ProductNote = ''.join(feature.qualifiers["note"]).replace(" [Streptomyces albus]", "").replace(" [Streptomyces albus J1074]", "")
				feature.qualifiers["note"] = ProductNote
			
			locusTag = ''.join(feature.qualifiers["locus_tag"])
			proteinID = "gnl|SalomonUMN|" + locusTag
			feature.qualifiers["protein_id"] = proteinID
			
			# feature.qualifiers.pop("EC_number", None)
			# feature.qualifiers.pop("function", None)
			
			product = ''.join(feature.qualifiers["product"])
			productStringLength = len(''.join(feature.qualifiers["product"]))
			if productStringLength > 100:
				print locusTag + "\t" + str(productStringLength) + "\t" + product
		
		if feature.type == "gene":
			feature.qualifiers.pop("function", None)
			# locusTag = str(feature.qualifiers["locus_tag"]).strip("[']").replace("Ssoud_0", "Ssoud_")
			# feature.qualifiers["locus_tag"] = locusTag
			
	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()