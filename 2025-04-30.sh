input="bigdata/2025-04-11/merged.tsv"
odir="results/2025-04-30"
o=$odir/README.md
mkdir -p $odir

annova(){
usage="$FUNCNAME <input.tsv> <odir/pref>"
if [ $# -lt 2 ];then echo "$usage";return;fi
input=$1
output=$2
mkdir -p ${output%/*}
Rscript -e '
    input="'$input'";
    o="'$output'";


    tt=read.table(input,header=T)
    p=tt[,grep(".numCs",colnames(tt))]/tt[,grep(".coverage",colnames(tt))]
    colnames(p)=gsub(".numCs",".percMeth",colnames(p))

    g=rep("Y",ncol(p)); g[grep("E",colnames(p))] = "E"; g[grep("W",colnames(p))] = "W"
    x= apply(p,1,function(x) sum(is.na(x)))
    tt=tt[x<3,]; p=p[x<3,]; ## filter out many n.a per row

    ## ANOVA
    pvals <- apply(p, 1, function(x) summary(aov(x ~ g,na.action = na.omit))[[1]]$`Pr(>F)`[1])
    pvals[is.na(pvals)]=1
    i= pvals < 0.05
    df=cbind(tt[,1:4],p,pvals)
    write.table(df,file=paste0(o,"_anova.tsv"),sep="\t",row.names=F,col.names=T,quote=F)

    library(ggplot2)
    library(reshape2)
    p_clean <- t(scale(t(p[i,])))  # or use: na.omit(p)
    p_clean$id <- 1:nrow(p_clean)
    meth_long <- melt(p_clean, id.vars = "id", variable.name = "Sample", value.name = "Methylation")
    ggplot(meth_long, aes(x = Sample, y = Methylation, fill = Sample)) +
      geom_violin(trim = FALSE, scale = "width") +
      theme_minimal() +
      labs(title = "Methylation Distribution per Sample", x = "Sample", y = "Methylation Level") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ## clustyering
    m=p[i,]
    m=t(scale(t(m)))
    h =Heatmap(m)
    row_dend <- row_dend(h)
    clusters <- cutree(as.hclust(row_dend), k = 5)
    png(file=paste0(o,"_heatmap.png"))
    top_anno = HeatmapAnnotation(foo = anno_boxplot(m, height = unit(2, "cm")))
    Heatmap(m,top_annotation= top_anno,row_split=clusters,show_row_names=F,cluster_columns=F,column_split=g)
    dev.off()

    for (k in unique(clusters)) {
      d <- df[i,][clusters == k, , drop = FALSE]
      write.table(d, file = paste0(o,"_cluster", k, ".tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
    }
'
}
pval=0.05
n=`tail -n+2 $odir/filt_anova.tsv | awk -v p=$pval '$(NF)<p' | wc -l`
N=`tail -n+2 $odir/filt_anova.tsv | wc -l `;
h=$odir/filt_heatmap.png;
echo '## Update : '`date`'
- Improved sensitivity in detecting '$n' / '$N' significant CpG sites (Pval <'$pval')
- Identified unbiased methylation patterns without applying prefiltering (previously required a 25% difference threshold)

## Results
- ANOVA results : [anova](/'$odir/filt_anova.tsv')
- Pval < '$pval' significant sites

| Heatmap | 
| :-: |
| ![hm](/'$h') |

| cluster_id | table |
| :-: | :-: |' > $o
for i in 1 2 3 5;do
    echo "| cluster_$i | [$i](/$odir/filt_cluster$i.tsv |" >> $o
done
git add -A
git commit -am wooutlier

exit
anno(){
i=filtered_anova.tsv
o=${i%.tsv}_anno
annotatePeaks.pl $i mm10 -go ${o}_go -annStats $o.stats > $o.tsv 
for i in 1 2 3 4 5;do
    continue;
    #i=merged_anova_pval005_cluster$i.tsv
    i=filtered_anova_cluster$i.tsv
    o=${i%.tsv}_anno
    annotatePeaks.pl $i mm10 -go ${o}_go -annStats $o.stats > $o.tsv 
done
}

bam2tsv(){
. ../../src/meth.sh
for f in *.bam;do
    n=${f%.dedu*}
    #bisbam2meth $f $n $n.tsv
    #chr1	4337982	4337982	+	11	8	3
done
}

cpg(){
    input=(
    2-year_1.tsv	2-year_3.tsv	E18pt5_1.tsv	E18pt5_3.tsv	Week4_1.tsv	Week4_3.tsv	
    2-year_2.tsv	2-year_4.tsv	E18pt5_2.tsv	E18pt5_4.tsv	Week4_2.tsv	Week4_4.tsv	
    )

    for f in ${input[@]};do
        n=${f%.tsv}
        #chr	start	end	strand	coverage	numCs	numTs
        tail -n+2 $f |  awk -v n=$n -v OFS="\t" '{print $1,$2,$2+1,n,$5,$6;}' |\
        intersectBed -a ../../data/cpgIslandExt.bed -b stdin -wa -wb |\
        sort -k1,1 -k2,3n | groupBy -g 1,2,3,8 -c 9,10, -o sum,sum
    done  |  awk '{print $1":"$2"-"$3"\t"$4"\t"($6/$5);}' | ca rc2mat -  > cpg.tsv
}

input=(
    E18pt5_1.tsv E18pt5_2.tsv E18pt5_3.tsv E18pt5_4.tsv
    Week4_1.tsv Week4_2.tsv Week4_3.tsv Week4_4.tsv
    2-year_1.tsv 2-year_2.tsv 2-year_3.tsv 2-year_4.tsv
)
. ../../src/util.sh;
#merge ${input[@]}  > merged.tsv

