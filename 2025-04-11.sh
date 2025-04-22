idir=bigdata/2025-04-11
odir=results/2025-04-11
mkdir -p $odir
cp $idir/*.png $odir/
cp $idir/filtered*.tsv $odir/
o=$odir/README.md

echo '

## Filtering Outliers 
- Outliers: E1, E2, W1, Y1, Y2

| before | after |
| :-: | :-: |
| ![sample_correlation](merged_anova_pval005_cor.png ) | ![filtered_correlation]( filtered_anova_cor.png ) |

## ANOVA and Clustering
- Used percent C methylation from filtered samples
- Applied a p-value threshold of 0.05 based on F-statistics for mean/variance differences (ANOVA)
- Grouped samples into 5 clusters based on methylation patterns

| Heatmap |
| :-: | 
|  ![filtered_heatmap]( filtered_anova_heatmap.png ) |
'> $o

echo '
## Annotation Per cluster

| cluster id | file |
| :-: | :-: |' >> $o
for f in $odir/filtered*anno.tsv;do
    n=` echo $f | grep -o "cluster\d" `
    echo "| $n | [$n.annotation](${f##*/}) |" >> $o
done

