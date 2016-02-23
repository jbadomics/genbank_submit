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
    description='corrects the format of Prokka-derived \inference strings in a GenBank file\nInspect your GenBank file and modify line 38 to match the accession numbers shown in \inference qualifiers\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file from which to extract hypothetical protein sequences')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

outputFile = open('output.gbk', 'w')

for sequenceRecord in SeqIO.parse(genbankFile, 'genbank'):

	for feature in sequenceRecord.features:
	
		if feature.type == "CDS":
		
			if len(feature.qualifiers["inference"]) >= 2:
				inferenceString = ''.join(feature.qualifiers["inference"][1])
				
				# Modify this search string as needed depending on your situation
				proteinAccession = re.search('\|YP_.*\|', inferenceString)
				
				if proteinAccession:
					newInferenceString = "similar to AA sequence:RefSeq:" + proteinAccession.group(0).replace("|", "")
					feature.qualifiers["inference"][1] = newInferenceString

	SeqIO.write(sequenceRecord, outputFile, "genbank")

outputFile.close()

genbankFile.close()
