#! /usr/bin/env pypy

# given a list of known index sequences and a FASTQ, count the number of times each index seems to be observed, with mismatch tolerance
# FASTQ sequences are truncated to the length of the known index sequences
# this will double-count reads that match more than one index!
# usage: collision_count.py index_table.tsv < index_reads.fastq
# index_table.tsv should be tab-delimited, two columns: first is name, second is sequence
# dumps frequency of known indices to STDOUT and most common reads (known or not) to STDERR

from __future__ import print_function
import sys, collections, distance

ALLOWED_MISMATCHES = 1 # bcl2fastq default
N_MOST_COMMON = 20 # number of most common reads to report

index_table = collections.OrderedDict()
index_counts = collections.Counter()
read_counts = collections.Counter()

for line in open(sys.argv[1]):
	name, sequence = line.rstrip().split('\t')
	index_table[name] = sequence

line_n = 0
for line in sys.stdin:
	line_n += 1
	if line_n % 4 == 1:
		assert(line.startswith('@'))
		index_read = line.rstrip()[line.rindex(':') + 1:]
		read_counts[index_read] += 1
		for name, known_sequence in index_table.items():
			index_counts[name] += (distance.hamming(known_sequence, index_read[:len(known_sequence)]) <= ALLOWED_MISMATCHES)

for name, sequence in index_table.items():
	print('%s\t%s\t%i' % (name, sequence, index_counts[name]))

for read_freq in read_counts.most_common(N_MOST_COMMON):
	print('%s\t%i' % read_freq, file = sys.stderr)

