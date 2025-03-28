#wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O data/mm10.chrom.sizes

for f in `ls -d /Volumes/T7/git/emseq/bigdata/*/`;do
    f=${f%/}
    s=${f##*/}
    echo $s
    ls -lah $f/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz 
done
#20250219_10_2-year_2_MY12882_S10_L007
#bismark/methylation_calls/methylation_coverage/sample_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz | awk '$5>10' | head
