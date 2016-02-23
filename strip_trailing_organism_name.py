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
    description='strips trailing organism name from protein products\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file to trim protein products')
parser.add_argument('[ORGANISM NAME]', help='full organism name to remove as it appears trailing some protein products;\ngenerally this is derived from a file of aa sequences used to annotate the genome in question\n*must be passed to the script in single quotes*')
parser.add_argument('[STRAIN]', help='full strain name to remove; generally this is derived from a file of aa sequences used to annotate the genome in question\n*must be passed to the script in single quotes*')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

organismName = sys.argv[2]
strain = sys.argv[3]

outputFile = open('output.gbk', 'w')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "CDS":
			
			fullString = " [" + organismName + strain + "]"
			fullStringNoSpace = "[" + organismName + strain + "]"
			organismOnlyString = " [" + organismName + "]"
			organismOnlyStringNoSpace = "[" + organismName + "]"
			product = ''.join(feature.qualifiers["product"]).replace(fullString, "").replace(fullStringNoSpace, "").replace(organismOnlyString, "").replace(organismOnlyStringNoSpace, "")
			feature.qualifiers["product"] = product
	
	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()