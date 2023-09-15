#! /usr/bin/env pypy3

# given a FASTQ file (stream) from Illumina including reads that failed the chastity filter,
# output only reads that passed the filter to STDOUT and report the counts of both to STDERR


import sys, re

# from Readfq by Heng Li, https://github.com/lh3/readfq
# yields reads as tuples of name, sequence, qualities (all strings)
def readfq(fp): # this is a generator function
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:], [], None # modified JWF to preserve entire read name
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield name, ''.join(seqs), None # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, seq, ''.join(seqs); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, seq, None # yield a fasta record instead
				break
# end of Readfq

def writefq (name, seq, qualities):
	return '@%s\n%s\n+\n%s\n' % (name, seq, qualities)

readname_regex = re.compile('^\S+ \d+:([YN]):\d+:[ACGT]+$')


n_pass = 0
n_fail = 0
for name, seq, qualities in readfq(sys.stdin):
	try:
		failed = readname_regex.match(name).groups()[0]
		if failed == 'N':
			n_pass += 1
			print(writefq(name, seq, qualities), end = '')
		elif failed == 'Y':
			n_fail += 1
		else:
			print('this should never happen')
			exit(1)
	except AttributeError:
		print('unrecognized read name format: ' + name, file = sys.stderr)
		exit(1)

print('passed: %i\nfailed: %i' % (n_pass, n_fail), file = sys.stderr)

