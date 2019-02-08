#! /bin/bash

# truncate a FASTQ file to the given length

if [ ! -n "$1" ]
then
	echo "usage: $(basename $0) LENGTH < INPUT_FASTQ > OUTPUT_FASTQ" >&2
	exit 1
fi

awk "{if (NR % 2 == 0) print substr(\$0, 1, $1); else print \$0}"

