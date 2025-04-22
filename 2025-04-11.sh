idir=bigdata/2025-04-11
odir=results/2025-04-11
mkdir -p $odir
cp $idir/*.png $odir/
cp $idir/filtered*.tsv $odir/
cp -r $idir/*_go $odir/

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

| cluster id | annotation file | GO |
| :-: | :-: | :-: |' >> $o
for f in $odir/filtered*anno.tsv;do
    go=https://github.com/hmgene/emseq/blob/main/${f%.tsv}_go/geneOntology.html
    go=https://raw.githack.com/hmgene/emseq/main/${f%.tsv}_go/geneOntology.html
    go="<a href=\"https://raw.githack.com/hmgene/emseq/main/${f%.tsv}/geneOntology.html\" target=\"_blank\"> View GO Results</a>"

    n=` echo $f | grep -o "cluster\d" `
    if [ -z "$n" ];then
        echo "| all | [all.annotation](${f##*/}) | $go |" >> $o
    else
        echo "| $n | [$n.annotation](${f##*/}) | $go |" >> $o
    fi
done

