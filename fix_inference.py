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
		
			if len(feature.qualifiers["inference"]) >= 2:
				inferenceString = ''.join(feature.qualifiers["inference"][1])
				proteinAccession = re.search('\|YP_.*\|', inferenceString)
				
				if proteinAccession:
					newInferenceString = "similar to AA sequence:RefSeq:" + proteinAccession.group(0).replace("|", "")
					print newInferenceString
					feature.qualifiers["inference"][1] = newInferenceString

	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()
