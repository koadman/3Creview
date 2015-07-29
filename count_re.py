#!/usr/bin/env python
# simple script to count up number of restriction sites in a sequence and its length
from __future__ import division
import fileinput
from Bio import SeqIO
import sys

if len(sys.argv) < 2:
	sys.exit("Usage: count_re.py <FastA sequence file>\n")

cutseq = 'CCGG'
for record in SeqIO.parse(open(sys.argv[1], "rU"), "fasta"):
	print record.id + "\t" + str(record.seq.count(cutseq) / len(record.seq)) + "\t" + str(len(record.seq))


