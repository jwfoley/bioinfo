#! /bin/bash

# given the path to an Illumina run folder and optionally a path to a sample sheet, run BCL Convert
# generates FASTQ files in the current directory and cleans up their names
# discards "Undetermined" reads by default unless "-u" is used
# passes any additional arguments to BCL Convert
# some useful BCL Convert arguments like "--no-lane-splitting" are better specified in the sample sheet ("NoLaneSplitting,true") to keep this script independent of library type and sequencer
# this version supports only single-end reads; fix this sometime

convert_cmd=bcl-convert
convert_args=( # hardcode additional settings here, e.g. thread numbers optimized for a specific computer
	'--output-directory .'
	'--force'
	'--fastq-gzip-compression-level 9'
	'--shared-thread-odirect-output true'
)
unwanted_name_regex='_S[0-9]+_R1_001' # fix this for R2

discard_undetermined=true
while getopts ":s:pu" opt; do
	case $opt in
		s)
			if [ ! -e "$OPTARG" ]; then echo "error: $OPTARG not found" >&2; exit 1; fi
			convert_args+=("--sample-sheet $OPTARG")
			;;
		p)
			convert_args+=('--bcl-sampleproject-subdirectories true')
			;;
		u)
			discard_undetermined=
			;;
	esac
done
shift "$((OPTIND-1))"
if [ $discard_undetermined ]; then
	convert_args+=('--bcl-only-matched-reads true')
fi


if [ ! -n "$1" ]; then
	echo "usage: $(basename $0) [-s sample_sheet] [-p] [-u] run_folder [...]
	-p: split output into subdirectories by Sample_Project
	-u: write Undetermined file
" >&2
	exit 1
fi
convert_args+=("--bcl-input-directory $1")
shift 1


set -euo pipefail


$convert_cmd ${convert_args[@]} "$@"

# clean up filenames
for fastq in $(ls *.fastq.gz */*.fastq.gz 2>/dev/null | grep -P $unwanted_name_regex); do
	mv "$fastq" $(sed -r "s/$unwanted_name_regex//" <<< "$fastq")
done

