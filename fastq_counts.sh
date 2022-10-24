#! /bin/bash

if [ ! -n "$1" ]
then
	echo "usage: $(basename $0) file1.fastq.gz file2.fastq.gz ..." >&2
	exit 1
fi

parallel 'echo -e {/.}\t$(( $(gzip -dc {} | wc -l) / 4 ))' ::: ${@[*]}

