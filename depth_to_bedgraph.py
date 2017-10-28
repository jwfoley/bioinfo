#! /usr/bin/env python3

# convert output of 'samtools depth' to UCSC bedGraph
# this is very simple and lets each position be a new entry

import sys

for line in sys.stdin:
	chr, pos, value = line.rstrip().split()
	if int(value) == 0: continue
	print('\t'.join([chr, str(int(pos) - 1), pos, value]))

