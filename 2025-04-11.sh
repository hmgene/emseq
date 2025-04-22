idir=bigdata/2025-04-11
odir=results/2025-04-11
mkdir -p $odir
cp $idir/*.png $odir/
cp $idir/filtered*.tsv $odir/
o=$odir/README.md
echo '
## Detect Outlier 

Filtered out E1, E2, W1, Y1, Y2

![sample_correlation](results/2025-04-11/merged_anova_pval005_cor.png )

## ANOVA and Clustering
- Used perent of C methylation of filtered samples
- 0.05 p-value threshold of F-statistics of mean/var difference
- 5 clusters given methylation patterns 

![filtered_correlation]( results/2025-04-11/filtered_anova_cor.png )

![filtered_heatmap]( results/2025-04-11/filtered_anova_heatmap.png )

## Annotation Per cluster
'> $o

for f in $odir/filtered*anno.tsv;do
    n=` echo $f | grep -o "cluster\d" `
    echo "- [$n]($f)" >> $o
done

