# MeersEtAl_MulTI-Tag
This repository contains code associated with the MulTI-Tag manuscript (Meers et al. 2021 bioRxiv):

Meers MP, Janssens DH, Henikoff S. Multifactorial chromatin regulatory landscapes at single cell resolution. _bioRxiv_ 2021.07.08.451691; DOI: [doi.org/10.1101/2021.07.08.451691](https://doi.org/10.1101/2021.07.08.451691)

This code is sufficient to generate processed single cell matrices underlying the UMAP and heatmap plots presented.

## Brief summary

1. **Align reads to hg19 genome build.** Fastq files are appended with barcode identity in the QNAME field for downstream processing into cell- and target-specific datasets.

2. **Merge aligned SAM files.** For downstream processing, SAM files representing H1 and K562 cell data from the same target should be merged using samtools merge.

3. **Generate CellRanger bed files.** We format data into CellRanger-style bed files of the following column structure:
    1) chr
    2) start
    3) end
    4) cell barcode
    5) number of duplicates

4. **Calculate unique fragments per cell.** For each SAM file (including merged SAM files), we generate two-column files representing unique fragments per cell and barcode assigned to that cell, sorted by descending number of unique fragments.

5. **Filter cells based on unique fragments.** For this analysis, we consider only cells that meet all of the following criteria:
    1) \> 500 unique H3K27me3 fragments
    2) \> 200 unique H3K4me2 fragments
    3) \> 200 unique H3K36me3 fragments

6. **Call peaks from aggregated data.** We use SEACR v1.4 with the following conditions: -n norm -m stringent -e 5

7. **Map single cell fragments onto peaks.** We use bedtools intersect -wao to quantify fragment overlap counts for each peak in each cell.

8. **Perform dimensionality reduction and plot data in UMAP form.** Cell-by-peak matrices are filtered by the number of cells reporting a fragment overlap, transformed by term frequency-inverse document frequency (TF-IDF) and log, subjected to Singular Value Decomposition (SVD) and most variable feature selection, and plotted in UMAP two-dimensional space.
9. 
