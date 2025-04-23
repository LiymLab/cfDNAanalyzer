`encode_unified_exclusion_list.bed`  
This is the encode unified GRCh38 exclusion list (https://www.encodeproject.org/files/ENCFF356LFX/)

`standard.chrom.sizes`  
Chromosome sizes for hg38

`alternative_haplotypes.bed`  
`centromeres.bed`  
`fix_patches.bed`  
`gaps.bed` 
Downloaded from UCSC table browser (https://genome.ucsc.edu/cgi-bin/hgTables)


`k100_minus_exclusion_lists.mappable_regions.bed`  
A bed file of all regions to be used for GC correction.  
Contains all intervals with perfect mappability (mappability=1) for 100bp reads obtained from the UCSC genome browser (https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw)
Additionally, all regions in `encode_unified_exclusion_list.bed` `alternative_haplotypes.bed` `centromeres.bed` `fix_patches.bed` and `gaps.bed` have been removed
