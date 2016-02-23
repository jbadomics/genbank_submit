#!/usr/bin/env python

from Bio import SeqIO
import sys
import re
import argparse
from argparse import RawTextHelpFormatter

# silence Biopython warnings
import warnings
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# script help and usage
parser=argparse.ArgumentParser(
    description='corrects typos and other user-specified grammatical errors\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file from which to extract hypothetical protein sequences')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

outputFile = open('output.gbk', 'w')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "CDS":
			
			#specify corrections to make here
			if "note" in feature.qualifiers:
				ProductNote = ''.join(feature.qualifiers["note"]).replace("domain containing", "domain-containing")
				feature.qualifiers["note"] = ProductNote
			
			Product = ''.join(feature.qualifiers["product"]).replace("domain containing", "domain-containing")
			feature.qualifiers["product"] = Product
			# feature.qualifiers.pop("EC_number", None)
			# feature.qualifiers.pop("function", None)
		
		if feature.type == "gene":
			feature.qualifiers.pop("function", None)
			# locusTag = str(feature.qualifiers["locus_tag"]).strip("[']").replace("Ssoud_0", "Ssoud_")
			# feature.qualifiers["locus_tag"] = locusTag
			
	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()