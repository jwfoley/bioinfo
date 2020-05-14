#! /usr/bin/env pypy3

# split a FASTA file containing multiple reference sequences into one file per sequence

import sys

out_file = None

for line in sys.stdin:
	if line.startswith('>'):
		out_file = open(line[1:].strip() + '.fa', 'w')
	out_file.write(line)

