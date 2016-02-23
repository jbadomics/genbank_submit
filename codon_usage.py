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
parser.add_argument('[CODON]', help='three-letter DNA codon of interest (e.g. ATG, not AUG)')
args=parser.parse_args()

codon = sys.argv[2]

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
					if 'ATA' in codonList:
						proteins_list.append(locusTag)
						outputString = '%s\t%i' % (locusTag, codonList.count(codon))
						outputFile.write(outputString + "\n")
print '%i proteins contain at least one %s codon' % (len(proteins_list), codon)
