# Integration with RNA-seq and ATAC-seq


## Data

| Type | File |
| :-: | :-: |
| RNA-seq LogFold | Beaudin.txt [rna_fl](../../data/from_pub/Beaudin.txt) |
| RNA-seq LogFold | Tan.txt [rna_bm](../../data/from_pub/Tan.txt) |
| ATAC-seq Fold | BM_HSC_on_FL_HSC_ATAC.xls [atac_bm](../../data/from_pub/BM_HSC_on_FL_HSC_ATAC.xls) |
| ATAC-seq Fold | FL_HSC_on_BM_HSC_ATAC.xls [atac_fl](../../data/from_pub/FL_HSC_on_BM_HSC_ATAC.xls) |
| EM-seq Methyl-Percent | filtered.3x.10bp.anova.anno.trend.tsv.gz [emseq_3x10bp](../2025-06-03-3x-depth/filtered.3x.10bp.anova.anno.trend.tsv.gz) |

## EM-seq data for IGV and GenomeBrowser 

| Type | File |
| :-: | :-: |
| bed Peak | Higher Methylation at W compared to E : [EtoW_up.bed](EtoW_up.bed)  |
| bed Peak | Lower Methylation at W compared to E  : [EtoW_dn.bed](EtoW_up.bed) |
| bed Peak | Higher Methylation at Y compared to W : [WtoY_up.bed](WtoY_up.bed) |
| bed Peak | Lower Methylation at Y compared to W : [WtoY_dn.bed](WtoY_up.bed) |
| bedGraph Averaged Percent | [mean_E.bedGraph.gz](mean_E.bedGraph.gz) |
| bedGraph Averaged Percent | [mean_W.bedGraph.gz](mean_W.bedGraph.gz) |
| bedGraph Averaged Percent | [mean_Y.bedGraph.gz](mean_Y.bedGraph.gz) | 
| bedGraph Averaged Percent | [meth_rna_atac.csv](meth_rna_atac.csv) |


## Heatmap Figures
- Significant genes were selected from ATAC-seq and RNA-seq datasets using a p-value threshold of 0.001.
- logFold and Fold values were rescaled column-wise.
- Heatmaps were generated using ComplexHeatmap with clustering applied.
- code :  [link](../../2025-07-14-public-data.sh)


| Ordering Criteria | Heatmap |
|-------------------|---------|
| By Methylation    | ![order_by_meth](meth_vs_rnaatac.heatmap.pdf) |
| By RNA            | ![order_by_rna](rnaatac_vs_meth.heatmap.pdf) |
| By Combined Signal| ![mix](meth_rna_atac_mix.heatmap.pdf) |



