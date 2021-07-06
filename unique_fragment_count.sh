#!/usr/bin/bash

## Step 4: Calculate unique fragments per cell

## NOTE: fastq files have been modified to contain barcodes delimited by "_" in the QNAME field
## in the following order: primer i7, primer i5, revcomp adapter i7, adapter i5

set -ue

if [ $# -lt 3 ]
then
	echo "
		Need at least three inputs: input bed file, output path, and output prefix!
	"
	exit 1
fi

i=`basename $1`
output="${2}/${3}.cellcount"

samtools view -bS -F 12 $1 | bedtools bamtobed  -bedpe -i - | awk '{split($7,s,"_"); print s[2]s[3]s[4]s[5]"\t"$1"\t"$2"\t"$6}' | sort | uniq | cut -f 1 | sort | uniq -c | sort -k1,1nr > $output

