## mapping summary

input="
Y1 bigdata/methylseq/2-year_1/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y2 bigdata/methylseq/2-year_2/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y3 bigdata/methylseq/2-year_3/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y4 bigdata/methylseq/2-year_4/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E1 bigdata/methylseq/E18pt5_1/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E2 bigdata/methylseq/E18pt5_2/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E3 bigdata/methylseq/E18pt5_3/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E4 bigdata/methylseq/E18pt5_4/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W1 bigdata/methylseq/Week4_1/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W2 bigdata/methylseq/Week4_2/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W3 bigdata/methylseq/Week4_3/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W4 bigdata/methylseq/Week4_4/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
"  
echo "sample group
E1 E
E2 E
E3 E
E4 E
W1 W
W2 W
W3 W
W4 W
Y1 Y
Y2 Y
Y3 Y
Y4 Y" > meta.txt

dn(){
    scp "hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/filtered.*.*bp.anova.tsv.gz" bigdata/
}

odir="results/2025-06-03-3x-depth"; mkdir -p $odir
out=$odir/README.md

cat bigdata/methylseq/*/out/bismark/summary/bismark_summary_report.txt | head -n 1 > $odir/summary.tsv
for f in bigdata/methylseq/*/out/bismark/summary/bismark_summary_report.txt;do
    s=`echo $f| cut -d"/" -f 3`;
    tail -n+2 $f | cut -f 2- | awk -v OFS="\t" -v s=$s '{ print s,$0;}'  >> $odir/summary.tsv
done
echo '## summary 
 [sum]('summary.tsv')

### BedGraph Files
'> $out

for x in 3x 5x;do
for b in 1bp 10bp;do
for i in E W Y;do
    o=bg/mean_${i}.$x.$b.bedGraph.gz;mkdir -p $odir/bg 
    #gunzip -dc bigdata/filtered.$x.$b.anova.tsv.gz | hm cutn - chrom,start,end,mean_${i} | tail -n+2 | awk -vOFS="\t" '{print $1,int($2),int($3),$4;}'  | gzip -c > $odir/$o
    echo " - [${o##*/}]($o) " >> $out
done
done
done

echo " 
3x10b bigdata/filtered.3x.10bp.anova.tsv.gz
5x10b bigdata/filtered.5x.10bp.anova.tsv.gz
3x1b bigdata/filtered.3x.1bp.anova.tsv.gz
5x1b bigdata/filtered.5x.1bp.anova.tsv.gz
" #| hm bismark-merge-anova -  $odir/multi_res_cor

echo '
## Multi resolution Correlation
 ![corr](multi_res_cor_heatmap.png)

## Trend Analysis

### Average Methylation %

| Group         | Sites | Mean E (%) | Mean W (%) | Mean Y (%) |
|---------------|-------|------------|------------|------------|
| 3prim         | 1405  | 63.68      | 62.18      | 57.69      |
| 5prim         | 238   | 54.65      | 55.01      | 48.14      |
| Intergenic    | 57719 | 66.99      | 64.80      | 65.03      |
| TTS           | 1813  | 63.25      | 61.61      | 59.46      |
| Exon          | 3712  | 55.99      | 62.62      | 65.01      |
| Intron        | 50188 | 66.36      | 62.11      | 56.78      |
| Non-coding    | 944   | 48.39      | 52.67      | 54.03      |
| Promoter-TSS  | 3910  | 53.80      | 54.20      | 53.03      |

###  E->W->Y Trends (E->W: up/dn, W->Y : up/dn)

| Group         | dn_dn | dn_nc | dn_up | nc_dn | nc_up | up_dn | up_nc | up_up |
|---------------|-------|-------|-------|-------|-------|-------|-------|-------|
| 3prim         | 373   | 0     | 356   | 2     | 1     | 379   | 4     | 290   |
| 5prim         | 66    | 1     | 52    | 0     | 0     | 68    | 1     | 50    |
| exon          | 659   | 1     | 896   | 2     | 1     | 1027  | 12    | 1114  |
| Intergenic    | 11706 | 40    | 18309 | 99    | 41    | 16211 | 211   | 11102 |
| intron        | 14279 | 49    | 13243 | 95    | 25    | 13484 | 145   | 8868  |
| non-coding    | 156   | 1     | 159   | 2     | 2     | 183   | 2     | 439   |
| promoter-TSS  | 863   | 3     | 1072  | 2     | 2     | 1140  | 13    | 815   |
| TTS           | 458   | 1     | 464   | 2     | 2     | 473   | 6     | 407   |

 [anova_annotation_trend.table]( filtered.3x.10bp.anova.anno.trend.tsv.gz )


![ trend ](vlnplot_methylation_trends.png)

![ trend_per_type]( vlnplot_methylation_trends_per_type.png)

![ go-term ](go-progressterm-per-trend.png)

| Trend    | Theme                                               | Stage                 | Key Systems                         |
|----------|-----------------------------------------------------|-----------------------|-------------------------------------|
| Up–Up    | early morphogenetic programs                        | Early → Late silencing | Neural crest, mesenchyme, skeletal  |
| Up–Down  | Temporary activation followed by repression         | Mid → Late            | Lung, brain, vasculature, synapse   |
| Down–Up  | Later activation for functional maturation          | Late                  | Synaptic, cardiac, skeletal         |



' >> $out


## run this on HPC
##echo "$input"| grep -v "^$" | hm bismark-merge-cov - 5 | gzip -c > bigdata/merged.x5.cov.gz 
##scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/merged.x5.cov.gz bigdata/
##scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/merged.3x.cov.gz bigdata/
##hm bismark-anova bigdata/filtered.3x.cov.gz meta.txt 0.05 > bigdata/filtered.cov.anova.tsv
exit

#hm bismark-cov-bin bigdata/filtered.cov.gz 1 | hm bismark-anova -  meta.txt 0.001 > bigdata/filtered.cov.${bz}bp.anova.p-3.tsv
#hm bismark-cov-groupfilter bigdata/merged.5x.cov.gz meta.txt 2 5 | gzip -c > bigdata/filtered.5x.cov.gz 
#hm bismark-cov-groupfilter bigdata/merged.5x.cov.gz meta.txt 2 5 | gzip -c > bigdata/filtered.3x.cov.gz 

exit
#hm bismark-cov-bin bigdata/filtered.cov.gz 10 | hm bismark-anova - meta.txt 0.05 > bigdata/filtered.cov.10b.anova.tsv
#for i in E W Y;do
#    hm cutn bigdata/filtered.cov.10b.anova.tsv chrom,start,end,mean_${i} | tail -n+2 | awk -vOFS="\t" '{print $1,int($2),int($3),$4;}'  > results/mean_${i}.bedGraph
#done
#
#target=( 
#Lin28b
#Hmga2
#Hic2
#Igf2bp2
#Igf2bp3
#)
#

