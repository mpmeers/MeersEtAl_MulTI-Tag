#!/usr/bin/bash

## Step 2: Merge SAM files from H1 and K562 cells for each target ##

set -ue

prefix1=`echo $1 | head -c -3`
prefix2=`echo $2 | head -c -3`

output="${prefix1}.${prefix2}.sam"

samtools merge $output $1 $2
