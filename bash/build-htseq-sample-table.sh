#!/bin/bash
htseqDir=$1
condition=$2

echo -e "sampleName\tfileName\tcondition"
for filename in `ls $htseqDir`; do
	echo -e "$filename\t$htseqDir$filename\t$condition"
done

#for the record, old way
# ls $htseqDir |  awk -F "." 'BEGIN { print "sampleName fileName condition"}
# {print $1" "$htseqDir"/"$0" "$condition}
# END {}' | sed "s/ /\\t/g"