#! /usr/bin/env python3

# given a GTF file of transcript annotations, return a BED file of introns
# usage: get_introns.py < transcripts.gtf > introns.bed
# warning: output is not (quite) sorted

import sys, genomerator

genes = genomerator.GtfStream(sys.stdin)
for gene in genes:
	for transcript in gene.data['transcript']:
		if not 'exon' in transcript.data or len(transcript.data['exon']) == 1: continue
		sorted_exons = sorted(transcript.data['exon'])
		for i in range(len(sorted_exons) - 1):
			intron = genomerator.GenomeFeature(
				reference_id =  transcript.reference_id,
				left_pos =      sorted_exons[i].right_pos + 1,
				right_pos =     sorted_exons[i + 1].left_pos - 1,
				is_reverse =    transcript.is_reverse
			)
			print(intron.bed(reference_names = genes.references, name= transcript.data['transcript_id']))

