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

for i in {1..5};do
#TermID	Term	Enrichment	logP	Genes in Term	Target Genes in Term	Fraction of Targets in Term	Total Target Genes	Total Genes	Entrez Gene IDs	Gene Symbols
#GO:0061646	positive regulation of glutamate neurotransmitter secretion in response to membrane depolarization	0.000343737107427677	-7.97563341535179	2	1	0.2	D
    cut -f 2,4,11 results/2025-04-11/*cluster$i*_go/biological_process.txt |\
    perl -ne 'chomp; my ($t,$p,$g)=split /\t/,$_;
    if(exp($p) < 0.005){ 
            $t=~s/regulation of/reg.of/g;
            $t=~s/positive reg.of/p.reg.of/g;
            $t=~s/negative reg.of/n.reg.of/g;

            my @tt=split/\s+/,$t;
            my @gg=split/,/,$g;
            map { 
                #print join("\t", join("_","cluster'$i'",@tt[0..2]), $_,1),"\n";
                print join("\t", join("_",,@tt[0..3]),"cluster'$i'", 1),"\n";
            } @gg;
    }
    ' 
done | ca rc2mat -  | tee $odir/go_bp.txt | Rscript -e 'tt=read.table("stdin",header=T,sep="\t")
library(ComplexHeatmap);
m=as.matrix(tt[,-1]);row.names(m)=tt[,1];
png("'$odir/go_bp_heatmap.png'",height=1200)
Heatmap(m)
dev.off();
'
echo '
## Functional Annotation Per Cluster
Summary
| GO process Summary Heatmap |
| :-: |
| ![goh]("'$odir/go_bp_heatmap.png'") |

Details
| cluster id | annotation file | GO |
| :-: | :-: | :-: |' >> $o

for f in $odir/filtered*anno.tsv;do
    go=https://github.com/hmgene/emseq/blob/main/${f%.tsv}_go/geneOntology.html
    go=https://raw.githack.com/hmgene/emseq/main/${f%.tsv}_go/geneOntology.html
    n=` echo $f | grep -o "cluster\d" `
    if [ -z "$n" ];then
        echo "| all | [all.annotation](${f##*/}) | [go]( $go ) |" >> $o
    else
        echo "| $n | [$n.annotation](${f##*/}) | [go]( $go ) |" >> $o
    fi
done

