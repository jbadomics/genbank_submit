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
    description='prints ncRNA features lacking required ncRNA_class descriptor for manual inspection\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file to search for problematic ncRNAs')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "ncRNA":
		
			if feature.qualifiers["ncRNA_class"] == None:
				locusTag = ''.join(feature.qualifiers["locus_tag"])
				print locusTag

genbankFile.close()
