input=(
results/E18pt5_vs_2-year_diff_25p_05q_hyper.tsv
results/E18pt5_vs_2-year_diff_25p_05q_hypo.tsv
results/E18pt5_vs_Week4_diff_25p_05q_hyper.tsv
results/E18pt5_vs_Week4_diff_25p_05q_hypo.tsv
results/Week4_vs_2-year_diff_25p_05q_hyper.tsv
results/Week4_vs_2-year_diff_25p_05q_hypo.tsv
)
fn(){
    #annotatePeaks.pl $1 mm10 > $1.tmp
    if [ -f $1.tmp ];then return;fi
    mv $2 $1.tmp
    R --no-save -e 'tt=read.table("'$1'",header=T)
    colnames(tt)=c("Chr","Start","End", colnames(tt)[4:ncol(tt)])
    tt1=read.csv("'$1.tmp'",header=T,sep="\\t");
    tt1$Start = tt1$Start -1;
    tt=merge(tt,tt1)
    write.table(file="'$2'",tt,col.names=T,row.name=F,quote=F,sep="\\t")
    '
    rm $1.tmp
};export -f fn;
#parallel fn {} {.}_anno.tsv ::: ${input[@]};


for f in ${input[@]};do
    s=${f#*/};s=${s%_diff*}
    h="hyper";
    if [ ${f/hypo/} != $f ];then
        h="hypo";
    fi
    tail -n+2 $f | awk -v OFS="\t" -v s=$s -v h=$h '{print $1":"$2"-"($2+1),s"."h,1}'
done | ca rc2mat -  > results/merged_diff.tsv
R --no-save -e  '
library(ComplexHeatmap)

tt=read.table("results/merged_diff.tsv",header=T)
m=as.matrix(tt[,-1])
row.names(m)= tt[,1]
n=m; n[n>0] =1

pdf("results/merged_diff.heatmap.pdf",height=14)
Heatmap(n, row_names_gp = gpar(fontsize = 6))
dev.off();
'

