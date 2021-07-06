#!/usr/bin/bash

## Step 7: Map unique single cell fragments onto called peaks

set -ue

if [ $# -lt 6 ]
then
	echo "
		Need five inputs in this order:
		
		1. H3K27me3 SEACR peak calls from Step 6
		2. H3K4me2 SEACR peak calls from Step 6
		3. H3K36me3 SEACR peak calls from Step 6
		4. H3K27me3 filtered CellRanger bed file from Step 5
		5. H3K4me2 filtered CellRanger bed file from Step 5
		6. H3K36me3 filtered CellRanger bed file from Step 5
		
	"
	exit 1
fi

module load bedtools

k27_out="${1}.singlecellmap.bed"
k4_out="${2}.singlecellmap.bed"
k36_out="${3}.singlecellmap.bed"

sbatch -n 1 --wrap="bedtools intersect -wao -a $1 -b $4 | sort -k1,1 -k2,2n -k3,3n -k10,10 | bedtools groupby -g 1,2,3,10,12 -c 3, -o count > $k27_out"
sbatch -n 1 --wrap="bedtools intersect -wao -a $2 -b $5 | sort -k1,1 -k2,2n -k3,3n -k10,10 | bedtools groupby -g 1,2,3,10,12 -c 3, -o count > $k4_out"
sbatch -n 1 --wrap="bedtools intersect -wao -a $3 -b $6 | sort -k1,1 -k2,2n -k3,3n -k10,10 | bedtools groupby -g 1,2,3,10,12 -c 3, -o count > $k36_out"
