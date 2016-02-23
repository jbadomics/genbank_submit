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
    description='adds NCBI protein IDs to a GenBank record\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file to add protein IDs to')
parser.add_argument('[PREFIX]', help='string to use for assigning protein IDs, which will have the form gnl|SmithUMN|locus_tag, where\nprefix = PILastnameUNIV')
args=parser.parse_args()

genbankFile = open(sys.argv[1], 'r')

proteinPrefix = sys.argv[2]

outputFile = open('output.gbk', 'w')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "CDS":

			locusTag = ''.join(feature.qualifiers["locus_tag"])
			proteinID = "gnl|" + proteinPrefix + "|" + locusTag
			feature.qualifiers["protein_id"] = proteinID
			
	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()
