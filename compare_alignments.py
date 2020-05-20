#! /usr/bin/env python3

import sys, argparse, pysam

def compare_read_positions (alignment1, alignment2):
	'''
	compare the positions of two alignments
	output them to the appropriate files if using those
	return the comparison result
	'''
	unmapped1 = alignment1 is None or alignment1.is_unmapped
	unmapped2 = alignment2 is None or alignment2.is_unmapped
	
	if unmapped1 and unmapped2:
		if outputs['concordant'] is not None: outputs['concordant'].write(alignment1)
		return 'concordant'
	elif not unmapped1 and unmapped2:
		if outputs['file1_only'] is not None: outputs['file1_only'].write(alignment1)
		return 'file1_only'
	elif unmapped1 and not unmapped2:
		if outputs['file2_only'] is not None: outputs['file2_only'].write(alignment2)
		return 'file2_only'
	elif (
		alignment1.is_reverse == alignment2.is_reverse and
		alignment1.reference_name == alignment2.reference_name and
		alignment1.reference_start - alignment1.query_alignment_start == alignment2.reference_start - alignment2.query_alignment_start # ignore soft-clipping differences
	):
		if outputs['concordant'] is not None: outputs['concordant'].write(alignment1)
		return 'concordant'
	else:
		if outputs['discordant'] is not None:
			outputs['discordant'].write(alignment1)
			outputs['discordant'].write(alignment2)
		return 'discordant'


parser = argparse.ArgumentParser(description = 'Given two unsorted SAM/BAM files generated from the same original FASTQ data, detect which reads were aligned at the same position in both SAM files, which were aligned in one but not the other, and which were aligned in different places. Input files should preferably include unaligned reads and contain the reads in roughly the same order in both files (probably the order of the original FASTQ file) but buffering allows some deviation at the expense of memory. Output file for concordant reads will only contain the alignments from file1. Output file for discordant reads will contain both alignments of each read.')
parser.add_argument('-c', '--concordant', metavar = 'BAMFILE', help = 'output for concordantly aligned reads')
parser.add_argument('-1', '--file1-only', metavar = 'BAMFILE', help = 'output for reads present only in file 1')
parser.add_argument('-2', '--file2-only', metavar = 'BAMFILE', help = 'output for reads present only in file 2')
parser.add_argument('-d', '--discordant', metavar = 'BAMFILE', help = 'output for discordantly aligned reads')
parser.add_argument('file1', metavar = 'BAMFILE')
parser.add_argument('file2', metavar = 'BAMFILE')
args = parser.parse_args()


input1 = pysam.Samfile(args.file1)
input2 = pysam.Samfile(args.file2)

outputs = {
	'concordant': None if args.concordant is None else pysam.Samfile(args.concordant, 'wb', header = input1.header),
	'file1_only': None if args.file1_only is None else pysam.Samfile(args.file1_only, 'wb', header = input1.header),
	'file2_only': None if args.file2_only is None else pysam.Samfile(args.file2_only, 'wb', header = input2.header),
	'discordant': None if args.discordant is None else pysam.Samfile(args.discordant, 'wb', header = input1.header)
}

# verify same references
if not input1.references == input2.references: print('warning: input files do not have matching reference names', file = sys.stderr)

# dicts of read positions, by name, from each file
# read is removed from dict when detected in the other file
# so these should stay small if the files are not in very different orders
buffer1 = {}
buffer2 = {}

counts = {
	'concordant': 0,
	'file1_only': 0,
	'file2_only': 0,
	'discordant': 0
}

input1_done = input2_done = False
while not (input1_done and input2_done):
	# iterate through each file, switching back and forth except when finding reads that were already identified in the other file
	# if the files are in exactly the same order then this will alternate back and forth keeping only one read position in memory at any time
	
	try:
		while True:
			alignment1 = next(input1)
			if alignment1.is_secondary or alignment1.is_supplementary: continue
			if alignment1.query_name in buffer2:
				counts[compare_read_positions(alignment1, buffer2[alignment1.query_name])] += 1
				del buffer2[alignment1.query_name]
			else:
				buffer1[alignment1.query_name] = alignment1
				break
	except StopIteration:
		input1_done = True
	
	try:
		while True:
			alignment2 = next(input2)
			if alignment2.is_secondary or alignment2.is_supplementary: continue
			if alignment2.query_name in buffer1:
				counts[compare_read_positions(buffer1[alignment2.query_name], alignment2)] += 1
				del buffer1[alignment2.query_name]
			else:
				buffer2[alignment2.query_name] = alignment2
				break
	except StopIteration:
		input2_done = True

if not (len(buffer1) == len(buffer2) == 0):
	print('warning: %i reads in file1 not found in file2 and %i reads in file2 not found in file1\nthese are treated as unaligned in the file where they are missing' % (len(buffer1), len(buffer2)), file = sys.stderr)
	for alignment1 in buffer1.values():
		counts[compare_read_positions(alignment1, None)] +=1
	for alignment2 in buffer2.values():
		counts[compare_read_positions(None, alignment2)] +=1

print('TOTAL\t%i' % sum(counts.values()))
for category, count in counts.items(): print('%s\t%i' % (category, count))

