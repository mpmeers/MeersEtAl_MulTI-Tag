#!/usr/bin/bash

## Step 1: Map raw fastq files to genome build ##

set -ue

if [ $# -lt 3 ]
then
	echo "Need input file list, output suffix, and genome build path!"
	exit 
fi
list=`cat $1`

names=`basename -a $list | awk 'BEGIN{FS="_"}; {print $3}' | sort -u`

module load Bowtie2


prefix=`basename -a $list | head -n 1 | cut -c 1-6`  #hg19 only
#prefix=`basename -a $list | tail -n 1 | cut -c 1-8` #hg19-mm10 hybrid
suffR1=`basename -a $list | head -n 1 | tail -c 10`
suffR2=`basename -a $list | sed -n '2p' | tail -c 10`
path=`dirname $list | head -n 1`
echo -e "$prefix\t$suffR1\t$suffR2\t$path"

for i in $names
do
	R1="${path}/${prefix}${i}${suffR1}"	
	R2="${path}/${prefix}${i}${suffR2}"	
	sbatch -n 1 --wrap="bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -x $3 -1 $R1 -2 $R2 -S $i.$2.sam"
done
