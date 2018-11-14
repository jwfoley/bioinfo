#! /usr/bin/env python3

import collections, sys, pysam

n_files = len(sys.argv) - 1
if n_files < 1:
	print('usage: %s file1.bam file2.bam ...' % sys.argv[0], file = sys.stderr)
	sys.exit(1)
counts = collections.defaultdict(lambda: [0] * n_files)

# count sizes
for file_number in range(n_files):
	for alignment in pysam.AlignmentFile(sys.argv[file_number + 1]):
		if not (
			alignment.is_unmapped or
			alignment.is_qcfail or
			alignment.is_secondary or
		  alignment.is_supplementary
		):
			counts[alignment.query_alignment_length][file_number] += 1

# print header
print('\t'.join(['aligned length'] + sys.argv[1:]))
for length in range(max(counts.keys())):
	print('\t'.join(map(str, [length] + counts[length])))

