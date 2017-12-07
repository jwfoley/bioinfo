#! /usr/bin/env pypy

# count the number of times each index sequence is be observed
# if a table of known sequences is provided, count hits to each one, with 1 mismatch tolerance (bcl2fastq default)
# if a sample sheet is also provided, label known sequences with sample IDs
# FASTQ sequences are truncated to the length of the known index sequences
# this will double-count reads that match more than one index!

# usage: collision_count.py [index_table.tsv] [SampleSheet.csv] < index_reads.fastq

# index_table.tsv should be tab-delimited, two columns: first is name, second is sequence
# dumps most common read counts to STDOUR and counts of known indices to STDERR


from __future__ import print_function
import sys, collections

ALPHABET = 'ACGTN'
N_MOST_COMMON = 20 # number of most common reads to report


read_counts = collections.Counter()
index_table = collections.OrderedDict()
index_rainbow_table = {}
index_counts = collections.Counter()
sample_table = collections.defaultdict(str)


if len(sys.argv) >= 2:
	for line in open(sys.argv[1]):
		name, sequence = line.rstrip().split('\t')
		index_table[name] = sequence
		index_rainbow_table[name] = set([sequence])
		for i in range(len(sequence)):
			for base in ALPHABET:
				index_rainbow_table[name].add(sequence[:i] + base + sequence[i + 1:])

if len(sys.argv) >= 3:
	in_data = False
	index_column = None
	for line in open(sys.argv[2]):
		if in_data:
			parsed_line = line.rstrip().split(',')
			sample_table[parsed_line[index_column]] = parsed_line[0]
		elif line.startswith('Sample_ID'):
			parsed_line = line.rstrip().split(',')
			column = 0
			for field in parsed_line:
				if field == 'index':
					index_column = column
					break
				column += 1
			in_data = True


line_n = 0
for line in sys.stdin:
	line_n += 1
	if line_n % 4 == 1:
		assert(line.startswith('@'))
		index_read = line.rstrip()[line.rindex(':') + 1:]
		read_counts[index_read] += 1
		for name, known_sequence in index_table.items():
			index_counts[name] += index_read[:len(known_sequence)] in index_rainbow_table[name]


for read_freq in read_counts.most_common(N_MOST_COMMON):
	print('%s\t%i' % read_freq)

for index_name, sequence in index_table.items():
	print('%s\t%s\t%s\t%i' % (sample_table[sequence], index_name, sequence, index_counts[index_name]), file = sys.stderr)	

