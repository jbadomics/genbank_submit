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
    description='collects statistics of codon presence/absence and total count per locus tag',
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nhttp://github.com/jbadomics/genbank_submit\nFebruary 2016\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file from which to extract genome sequence')
parser.add_argument('[CODON]', help='three-letter DNA codon of interest (e.g. ATG, not AUG)')
args=parser.parse_args()

codon = sys.argv[2].upper().replace("U", "T")

if len(codon) % 3 != 0:
	print 'codon provided is not a multiple of 3'
	sys.exit()

outputFileName = sys.argv[1].replace(".gbk", "." + codon + ".txt")
proteins_list = []
with open(outputFileName, 'wb') as outputFile:
	with open(sys.argv[1], 'r') as genbankFile:
		for sequenceRecord in SeqIO.parse(genbankFile, "genbank"):
			genome = sequenceRecord.seq
			for feature in sequenceRecord.features:
				strand = int(feature.location.strand)
				if feature.type == 'CDS':
					if strand == 1:
						locusTag = ''.join(feature.qualifiers["locus_tag"])
						startCoord = int(feature.location.start.position)
						endCoord = int(feature.location.end.position)
						geneSequence = str(genome[startCoord:endCoord])
						codonList = [geneSequence[i:i+3] for i in range(0, len(geneSequence))]
					if strand == -1:
						locusTag = ''.join(feature.qualifiers["locus_tag"])
						startCoord = int(feature.location.end.position) #this is not a typo
						endCoord = int(feature.location.start.position) #this is not a typo
						geneSequence = str(genome[endCoord:startCoord].reverse_complement())
						codonList = [geneSequence[i:i+3] for i in range(0, len(geneSequence))]
					if codon in codonList:
						proteins_list.append(locusTag)
						outputString = '%s\t%i' % (locusTag, codonList.count(codon))
						outputFile.write(outputString + "\n")
if len(proteins_list) > 0:
	print '%i proteins contain at least one %s codon' % (len(proteins_list), codon)
	print 'results written to %s' % (outputFileName)
else:
	print '%s not found' % (codon)

