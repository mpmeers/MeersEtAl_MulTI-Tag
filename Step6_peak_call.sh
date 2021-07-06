#!/usr/bin/bash

## Step 6: Convert filtered CellRanger bed files to bedGraphs and call peaks using SEACR

## NOTE: bedtools and SEACR v1.4 must be in your path

set -ue

if [ $# -lt 4 ]
then
	echo "
		Need five inputs in this order:
		
		1. H3K27me3 CellRanger unique fragment count-filtered bed file
		2. H3K4me2 CellRanger unique fragment count-filtered bed file
		3. H3K36me3 CellRanger unique fragment count-filtered bed file
		4. Genome size file (for bedtools genomecov). Format: <chr>\t<chrlength>
		5. IgG control bedgraph file (for SEACR peak calling)
		
	"
	exit 1
fi

module load bedtools

sbatch -n 1 --wrap="bedtools genomecov -bg -i $1 -g $4 | bash SEACR_1.4.sh - $5 norm stringent $1 5"
sbatch -n 1 --wrap="bedtools genomecov -bg -i $2 -g $4 | bash SEACR_1.4.sh - $5 norm stringent $2 5"
sbatch -n 1 --wrap="bedtools genomecov -bg -i $3 -g $4 | bash SEACR_1.4.sh - $5 norm stringent $3 5"
