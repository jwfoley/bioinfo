# forward = R1 CT, R2 GA
# reverse = R1 GA, R2 CT
# for convenience, this script only looks at R1

echo -e 'library\tforward\treverse'
for bam in $@
do
	echo -e $(basename $bam) '\t' $(samtools view -f 0x40 $bam | grep -c "ZB:Z:CT$") '\t' $(samtools view -f 0x40 $bam | grep -c "ZB:Z:GA$")
done

