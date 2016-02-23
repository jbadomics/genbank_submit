#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
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
    description='extracts amino acid sequences of hypothetical proteins FASTA format from a (multi)-Genbank file.',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/tnseq\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file from which to extract genome sequence')
args=parser.parse_args()

outputFileName = sys.argv[1].replace(".gbk", ".hypotheticals.faa")

with open(outputFileName, 'wb') as outputFile:
	with open(sys.argv[1], 'r') as genbankFile:
		for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
			for feature in sequenceRecord.features:
				if feature.type == 'CDS':
					if 'ypothetical' in ''.join(feature.qualifiers["product"]):
						locusTag = ''.join(feature.qualifiers["locus_tag"])
						aaSeq = ''.join(feature.qualifiers["translation"])
						outputString = '>%s\n%s' % (locusTag, aaSeq)
						#print outputString
						#SeqIO.write(outputString, outputFile, "fasta")
						outputFile.write(outputString + "\n")
