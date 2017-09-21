#!/usr/bin/python3

# get ERCC hit counts from a list of BAM files and print a nice matrix
# usage: ercc_counts.py file1.bam file2.bam ...

import pysam, sys

bamfiles = sys.argv[1:]

counts = {}

for bamfile_index in range(0, len(bamfiles)):
	sam = pysam.Samfile(bamfiles[bamfile_index], 'rb')
	which_ercc = [tid for tid in range(0, sam.nreferences) if sam.getrname(tid).startswith('ERCC')]
	
	for ercc_tid in which_ercc:
		ercc_name = sam.getrname(ercc_tid)
		if not ercc_name in counts:
			counts[ercc_name] = [0] * len(bamfiles)
	
	for read in sam:
		if read.is_secondary or read.is_unmapped or read.is_read2:
			continue
		if read.reference_id in which_ercc:
			counts[sam.getrname(read.reference_id)][bamfile_index] += 1

print('\t'.join([''] + bamfiles))
for ercc_name in sorted(counts.keys()):
	print('\t'.join([ercc_name] + [str(x) for x in counts[ercc_name]]))


