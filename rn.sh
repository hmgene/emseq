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
"""
           group  Sites    mean_E     mean_W     mean_Y
          <char> <int>      <num>      <num>      <num>
            <NA>    14  0.2340333  0.2077105  0.5092456
          3prim  1405 63.6838120 62.1841517 57.6932052
          5prim  238 54.6511712 55.0127968 48.1433496
      Intergenic 57719 66.9876675 64.7996101 65.0313108
             TTS  1813 63.2491014 61.6091227 59.4552735
            exon  3712 55.9893814 62.6223668 65.0061654
          intron 50188 66.3645757 62.1072053 56.7808921
      non-coding   944 48.3941965 52.6688076 54.0259305
    promoter-TSS  3910 53.7976471 54.1951853 53.0326458
"""

###  E->W->Y Trends (E->W: up/dn, W->Y : up/dn)
"""
            group dn_dn dn_nc dn_up nc_dn nc_up up_dn up_nc up_up
            3prim   373     0   356     2     1   379     4   290
            5prim    66     1    52     0     0    68     1    50
             exon   659     1   896     2     1  1027    12  1114
       Intergenic 11706    40 18309    99    41 16211   211 11102
           intron 14279    49 13243    95    25 13484   145  8868
       non-coding   156     1   159     2     2   183     2   439
     promoter-TSS   863     3  1072     2     2  1140    13   815
              TTS   458     1   464     2     2   473     6   407
             <NA>     1     0     7     0     0     1     0     5
"""

 [anova_annotation_trend.table]( filtered.3x.10bp.anova.anno.trend.tsv.gz )
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

