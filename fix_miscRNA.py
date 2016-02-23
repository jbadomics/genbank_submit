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
    description='collects statistics of codon presence/absence and total count per locus tag',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file to correct ydaO-yuaA miscRNA features')
args=parser.parse_args()
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
