#! /usr/bin/env python3

# given a specific sequence (e.g. restriction enzyme's motif), make a BED file of all perfect matches in a reference genome (FASTA)
# allow target sequence or reference genome to contain ambiguous bases in IUPAC format
# this probably only works if the sequence is a palindrome (like a restriction enzyme's motif) and it won't find overlapping hits

import re, sys
from Bio import SeqIO

N_LIMIT = 1 # don't allow hits with more than this many Ns
IUPAC = { # each key is a base code and each value is *all the base codes that it can contain*
	'R': 'AG',
	'Y': 'CT',
	'S': 'GC',
	'W': 'AT',
	'K': 'GT',
	'M': 'AC',
	'B': 'CGT',
	'D': 'AGT',
	'H': 'ACT',
	'V': 'ACG',
	'N': 'ACGT'
}
def base_to_regex (base):
	try:
		these_bases = IUPAC[base]
		result = these_bases
		for code, bases in IUPAC.items():
			if set(these_bases).intersection(bases):
				result += code
	except KeyError:
		result = base
		for code, bases in IUPAC.items():
			if base in bases: result += code
	return('[' + ''.join(sorted(result)) + ']')

if len(sys.argv) < 2: sys.exit('usage: %s [target sequence] < [reference FASTA] > [BED file]' % sys.argv[0])

# convert target sequence into regex
target_re = re.compile(''.join(list(map(base_to_regex, sys.argv[1].upper()))))


for chrom in SeqIO.parse(sys.stdin, 'fasta'):
	for hit in target_re.finditer(str(chrom.seq.upper())):
		if hit.group().count('N') <= N_LIMIT:
			if hit.group().count('R') > 0: print('\t'.join((chrom.name, str(hit.start()), str(hit.end()))))

