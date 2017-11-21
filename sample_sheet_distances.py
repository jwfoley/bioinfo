#! /usr/bin/env python3

# given an Illumina sample sheet (as STDIN), return a distance matrix of the samples' barcode sequences
# currently only works on single-indexed libraries

import sys

def hamming_distance (seq1, seq2):
	assert len(seq1) == len(seq2)
	return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))

for line in sys.stdin:
	if line.rstrip() == 'Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description': break # header


indices = {}
for line in sys.stdin:
	parsed_line = line.rstrip().split(',')
	sample_id, index = parsed_line[0], parsed_line[5]
	assert sample_id not in indices
	indices[sample_id] = index

print('\t'.join([''] + list(indices.keys())))
for row_sample, row_sequence in indices.items():
	out_line = row_sample
	for col_sample, col_sequence in indices.items():
		out_line += '\t'
		if col_sample != row_sample:
			out_line += str(hamming_distance(row_sequence, col_sequence))
	print(out_line)

