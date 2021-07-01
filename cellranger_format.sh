#!/usr/bin/bash

## Step 2: Convert mapped SAM files to CellRanger-formatted fragment bed files:

## 1) chr
## 2) start
## 3) end
## 4) cell barcode
## 5) fragment duplicates

set -ue

if [ $# -lt 3 ]
then
	echo "
		Need at least three inputs: input bed file, output path, and output prefix!
	"
	exit 1
fi

output="${2}/${3}.sci.bed"

module load samtools
module load bedtools

samtools view -b -F 12 $1 | bedtools bamtobed -bedpe -i stdin | awk '$1==$4 && $6-$2 <= 1000 {split($7, s, "_"); print $1"\t"$2"\t"$6"\t"s[2]s[3]s[4]s[5]}' | sort -k4,4 -k1,1 -k2,2n -k3,3n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | sort -k1,1 -k2,2n -k3,3n > $output
