#! /usr/bin/env python3

# given a list of strings of the same length (e.g. index sequences), report a Hamming distance matrix
# remember you can terminate console STDIN with Ctrl-D
# writes distance matrix to STDOUT and summary of distances to STDERR

from scipy.spatial.distance import hamming
from collections import Counter
import sys

sequences = []
seq_length = None
for line in sys.stdin:
	seq = line.strip()
	if len(seq) == 0: continue
	
	# keep all indexes the same length
	if seq_length is None:
		seq_length = len(seq)
	elif len(seq) != seq_length:
		print('error: different length sequence %s' % seq, file = sys.stderr)
		exit(1)
	
	sequences += [seq]

print('\t'.join([''] + sequences))
distances = Counter()
for i in range(len(sequences)):
	output = [sequences[i]]
	for j in range(0, i):
		dist = int(seq_length * hamming(tuple(sequences[i]), tuple(sequences[j])))
		output += [str(dist)]
		distances[dist] += 1
	print('\t'.join(output))

print('distance\tfrequency', file = sys.stderr)
for dist, count in sorted(distances.items()): print('%i\t%i' % (dist, count), file = sys.stderr)

