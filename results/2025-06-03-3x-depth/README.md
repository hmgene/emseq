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


