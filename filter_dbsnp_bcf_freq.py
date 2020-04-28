#! /usr/bin/env python3

# given a BCF file of dbSNP variant calls in NCBI format, produce a BCF file with only variants above a minimum minor-allele frequency
# NCBI format shows results from multiple studies and a variant is kept if any of them reported a minor-allele frequency above the threshold
# where more than one alternate allele was given, only those that pass the frequency threshold are reported, so any given variant position from the input might report fewer alternate alleles in the output (but the FREQ information is unaltered, making it inconsistent!)
# input is automatically recognized from BCF or VCF
# output is in VCF format because that's what I needed when I wrote this

import sys, pysam

if len(sys.argv) != 2: sys.exit('usage: %s MIN_FREQ < IN.BCF > OUT.VCF' % sys.argv[0])
freq_threshold = float(sys.argv[1])

in_file = pysam.VariantFile('-')
out_file = pysam.VariantFile('-', 'w', header = in_file.header)
for variant in in_file:
	try:
		# NCBI format doesn't parse cleanly, so the first element is the source name plus the major-allele frequency, then each successive element is a minor-allele frequency possibly concatenated with the next source's name and major-allele frequency 
		passed = [False] * len(variant.alts) # whether each alternate allele passes the threshold (maybe one does and another doesn't)
		for freq_str, alt_index in zip(variant.info['FREQ'][1:], range(len(variant.alts))):
			try:
				if float(freq_str.partition('|')[0]) >= freq_threshold: passed[alt_index] = True
			except ValueError: # unable to recast as float: no frequency given
				continue
		if sum(passed) > 0: # any alternates passed
			variant.alts = ((alt for alt, alt_passed in zip(variant.alts, passed) if alt_passed)) # select only alternates that passed the threshold
			out_file.write(variant)
	except KeyError: # no FREQ so it automatically fails
		continue

