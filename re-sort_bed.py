#! /usr/bin/env python3

# given a BED file and a list of reference names, sort the BED file according to the provided order of reference names

import sys, genomerator

if len(sys.argv) != 2:
	print('usage: %s [list of reference names] < [input BED] > [output BED]' % sys.argv[0])
	exit(1)

references = genomerator.read_references(open(sys.argv[1]))

input_bed = genomerator.BedStream(sys.stdin, references = references)

try:
	for feature in sorted(input_bed): print(feature.bed(reference_names = references))
except KeyError:
	pass


