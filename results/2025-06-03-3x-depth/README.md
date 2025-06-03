## summary 
 [sum](summary.tsv)

### BedGraph Files

 - [mean_E.3x.1bp.bedGraph.gz](bg/mean_E.3x.1bp.bedGraph.gz) 
 - [mean_W.3x.1bp.bedGraph.gz](bg/mean_W.3x.1bp.bedGraph.gz) 
 - [mean_Y.3x.1bp.bedGraph.gz](bg/mean_Y.3x.1bp.bedGraph.gz) 
 - [mean_E.3x.10bp.bedGraph.gz](bg/mean_E.3x.10bp.bedGraph.gz) 
 - [mean_W.3x.10bp.bedGraph.gz](bg/mean_W.3x.10bp.bedGraph.gz) 
 - [mean_Y.3x.10bp.bedGraph.gz](bg/mean_Y.3x.10bp.bedGraph.gz) 
 - [mean_E.5x.1bp.bedGraph.gz](bg/mean_E.5x.1bp.bedGraph.gz) 
 - [mean_W.5x.1bp.bedGraph.gz](bg/mean_W.5x.1bp.bedGraph.gz) 
 - [mean_Y.5x.1bp.bedGraph.gz](bg/mean_Y.5x.1bp.bedGraph.gz) 
 - [mean_E.5x.10bp.bedGraph.gz](bg/mean_E.5x.10bp.bedGraph.gz) 
 - [mean_W.5x.10bp.bedGraph.gz](bg/mean_W.5x.10bp.bedGraph.gz) 
 - [mean_Y.5x.10bp.bedGraph.gz](bg/mean_Y.5x.10bp.bedGraph.gz) 

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

