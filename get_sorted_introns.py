#! /usr/bin/env python3

# given a sorted GTF file of transcript annotations, return a sorted BED file of introns
# usage: get_introns.py < transcripts.gtf > introns.bed

import sys, genomerator

gene_input = genomerator.GtfStream(sys.stdin, verify_sorted = True)
gene_overlapper = genomerator.FeatureOverlapper(gene_input)
for genes in gene_overlapper:
	intron_buffer = []
	for gene in genes:
		for transcript in gene.data['transcript']:
			if not 'exon' in transcript.data or len(transcript.data['exon']) == 1: continue
			sorted_exons = sorted(transcript.data['exon']) # in case they're in reverse or something
			for i in range(len(sorted_exons) - 1):
				intron_buffer.append(genomerator.GenomeFeature(
					reference_id =  transcript.reference_id,
					left_pos =      sorted_exons[i].right_pos + 1,
					right_pos =     sorted_exons[i + 1].left_pos - 1,
					is_reverse =    transcript.is_reverse,
					data =          transcript.data['transcript_id']
				))
	
	# re-sort the introns to be safe
	for intron in sorted(intron_buffer):
		print(intron.bed(reference_names = gene_input.references, name = intron.data))

