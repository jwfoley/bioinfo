#! /usr/bin/env python3

import sys, os, argparse, gffutils, pyfaidx

DB_FILE = '/tmp/transcripts_to_fasta.db'

parser = argparse.ArgumentParser(description = 'given a GFF/GTF file of transcript annotations and a FASTA file of genome sequence, produce a strand-correct FASTA file of transcript sequences')
parser.add_argument('fasta_file', action = 'store')
parser.add_argument('gtf_file', action = 'store')
parser.add_argument('output_file', type = argparse.FileType('w'), default = sys.stdout, nargs = '?')
parser.add_argument('-i', '--infer', action = 'store_true', help = 'infer transcripts (not listed separately in GTF')
args = parser.parse_args()

db = gffutils.create_db(
	args.gtf_file,
	DB_FILE,
	force = True,
	disable_infer_transcripts = not args.infer,
	disable_infer_genes = not args.infer
)

fasta = pyfaidx.Fasta(args.fasta_file)
for transcript in db.features_of_type('transcript'):
	transcript_id_list = transcript['transcript_id']
	assert(len(transcript_id_list) == 1)
	transcript_id = transcript_id_list[0]
	
	seq = ''
	for exon in db.children(
		transcript_id,
		featuretype = 'exon',
		order_by = 'start',
		reverse = (transcript.strand == '-')
	):
		assert exon.chrom == transcript.chrom
		seq += exon.sequence(fasta)
	
	args.output_file.write('>%s\n%s\n' % (transcript_id, seq))

os.remove(DB_FILE)

