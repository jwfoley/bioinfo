#! /bin/bash

tmp_dir=/tmp
bedgraph_to_bigwig=~/bin/bedGraphToBigWig
depth_to_bedgraph="$(dirname $0)/depth_to_bedgraph.py"

region_arg=''
while getopts ":r:" opt
do
	case $opt in
		r)
			region_arg="-R $OPTARG"
			;;
	esac
done
shift "$((OPTIND-1))"

if [ ! -n "$3" ]
then
	echo "usage: $(basename $0) [-r region] chrom.sizes output_name file1.bam file2.bam ..." >&2
	exit 1
fi

chrom_sizes=$(readlink -f $1)
output_name=$2
shift 2


set -euo pipefail

samtools merge $region_arg -u - "$@" | tee >(samtools view -uF 0x10 - | samtools depth - | $depth_to_bedgraph > $tmp_dir/tmp+.bedgraph ) | samtools view -uf 0x10 - | samtools depth - | sed 's/\(.*\t\)/\1-/' | $depth_to_bedgraph > $tmp_dir/tmp-.bedgraph # the sed command puts a minus sign before values on the minus strand


# bedGraphToBigWig requires seekable files for some reason
$bedgraph_to_bigwig $tmp_dir/tmp+.bedgraph $chrom_sizes $output_name+.bw
rm $tmp_dir/tmp+.bedgraph
$bedgraph_to_bigwig $tmp_dir/tmp-.bedgraph $chrom_sizes $output_name-.bw
rm $tmp_dir/tmp-.bedgraph
