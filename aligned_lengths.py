#! /usr/bin/env python3

import collections, sys, pysam, click

@click.command()
@click.option('-f', '--feature_counts', is_flag = True, help = 'only count feature- (gene-)aligned reads; second half of input files are the .featureCounts read-by-read summary files corresponding to the first half of input files')
@click.argument('files', nargs = -1)
def aligned_lengths(feature_counts, files):
	if feature_counts:
		if not len(files) % 2 == 0:
			raise RuntimeError('--feature_counts requires a list of BAM files and then a list of .featureCounts files in corresponding order')
		n_files = int(len(files) / 2)
		input_files = files[:n_files]
		hit_files = files[n_files:]
	else:
		n_files = len(files)
		input_files = files
	
	counts = collections.defaultdict(lambda: [0] * n_files)

	# count sizes
	for file_number in range(n_files):
		if feature_counts: hit_stream = open(hit_files[file_number])
		for alignment in pysam.AlignmentFile(input_files[file_number]):
			if feature_counts:
				hit_line = next(hit_stream).split('\t')
				if not hit_line[0] == alignment.query_name: raise RuntimeError('found %s in %s but %s in %s' % (alignment.query_name, input_files[file_number], hit_line[0], hit_files[file_number]))
			if not (
				alignment.is_unmapped or
				alignment.is_qcfail or
				alignment.is_secondary or
				alignment.is_supplementary
			):			
				if (not feature_counts) or hit_line[1] == 'Assigned':
					counts[alignment.query_alignment_length][file_number] += 1

	# print header
	print('\t'.join(['aligned length'] + list(input_files)))
	# print table
	for length in range(max(counts.keys())):
		print('\t'.join(map(str, [length] + counts[length])))

if __name__ == '__main__':
	aligned_lengths()
