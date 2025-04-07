wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O data/mm10.chrom.sizes

samtools view -bh  bigdata/bams/20250219_10_2-year_2_MY12882_S10_L006.deduplicated.sorted.bam chrX:1-10000000 > data/sample.bam 

#/Users/hyunminkim/miniforge3/envs/r4.4/lib/R/library/methylKit/extdata/cpgi.hg18.bed.txt
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz -O data/cpgIslandExt.txt.gz
gunzip -dc  data/cpgIslandExt.txt.gz  | awk -v OFS="\t" '{ print $2,$3,$4,$5"_"$6;}' > data/cpgIslandExt.bed

cp /Volumes/T7/git/emseq/bigdata/20250219_10_2-year_2_MY12882_S10_L006/bismark/methylation_calls/methylation_coverage/sample_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz bigdata/covs/


