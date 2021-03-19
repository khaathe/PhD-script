#!/bin/bash
dir="./Result/bowtie2/alignment/GRCh38_copy"

for sam in `find $dir -iname "*\.sam" `; do
	bam=`echo $sam|sed -e "s/.sam/.bam/g"`
	samtools view -b $sam -o $bam
done 