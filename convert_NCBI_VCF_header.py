#! /usr/bin/env python3

# convert the chromosome names in an NCBI-style VCF from RefSeq accession numbers to UCSC-style names
# use an NCBI assembly report file to look up the conversions
# if the report gives "na" for the UCSC-style name, keep the accession number
# can be used on an entire VCF but it's probably faster to just make a new header and use "bcftools reheader"

import sys

if len(sys.argv) != 2: sys.exit('%s assembly_report.txt < oldvcf > newvcf' % sys.argv[0])

conversion_table = {}
accession_col, name_col = None, None
for line in open(sys.argv[1]):
	if line.startswith('# Sequence-Name'):
		fields = line.rstrip().split('\t')
		accession_col = fields.index('RefSeq-Accn')
		name_col = fields.index('UCSC-style-name')
	elif accession_col is None or name_col is None:
		continue
	else:
		values = line.rstrip().split('\t')
		name = values[name_col]
		if name != 'na': conversion_table[values[accession_col]] = name

for line in sys.stdin:
	if line.startswith('##contig'):
		accession = line.rstrip()[line.find('ID=') + 3:-1]
		try:
			print('##contig=<ID=' + conversion_table[accession] + '>')
		except KeyError:
			print(line, end = '')
	else:
		print(line, end = '')

