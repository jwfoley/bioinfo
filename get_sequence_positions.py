#! /usr/bin/env python3

# given a specific sequence (e.g. restriction enzyme's motif), make a BED file of all perfect matches in a reference genome (FASTA)
# allow target sequence or reference genome to contain ambiguous bases in IUPAC format
# this will find overlapping hits
# if the target sequence isn't a palindrome, it will find hits on both strands and label the BED accordingly
# is there anything more useful to put in the BED?
# doesn't count hits with too many N in the reference sequence
# this is slow; use pypy


import regex, sys
from Bio.Seq import Seq
from Bio import SeqIO

N_LIMIT = 1 # don't allow hits with more than this many Ns in the reference genome
IUPAC = { # each key is a base code and each value is *all the base codes that it can contain*
	'A': 'A',
	'C': 'C',
	'G': 'G',
	'T': 'T',
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
	these_bases = IUPAC[base]
	result = ''
	for code, bases in IUPAC.items():
		if set(these_bases).intersection(bases):
			result += code
	return(''.join(sorted(result)))

if len(sys.argv) < 2: sys.exit('usage: %s [target sequence] < [reference FASTA] > [BED file]' % sys.argv[0])

# convert target sequence into regex
target_seq = Seq(sys.argv[1].upper())
target_seq_revcomp = target_seq.reverse_complement()
is_palindrome = target_seq_revcomp == target_seq
target_re = regex.compile(''.join(['[%s]' % bases for bases in list(map(base_to_regex, target_seq))]))
target_revcomp_re = regex.compile(''.join(['[%s]' % bases for bases in list(map(base_to_regex, target_seq_revcomp))]))

for chrom in SeqIO.parse(sys.stdin, 'fasta'):
	hits = []
	for hit in target_re.finditer(str(chrom.seq.upper()), overlapped = True):
		if hit.group().count('N') <= N_LIMIT:
			hits.append([chrom.name, hit.start(), hit.end(), '+']) # last value is strand
	if not is_palindrome:
		for hit in target_revcomp_re.finditer(str(chrom.seq.upper()), overlapped = True):
			if hit.group().count('N') <= N_LIMIT:
				hits.append([chrom.name, hit.start(), hit.end(), '-'])
		hits = sorted(hits, key = lambda x: x[1]) # try bisect.insort or bisect.insort_left
	
	for hit in hits:
		if is_palindrome:
			print('%s\t%i\t%i' % tuple(hit[:3]))
		else:
			print('%s\t%i\t%i\t.\t.\t%s' % tuple(hit))

