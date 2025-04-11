idir="../results/2025-04-04"
iter=( "E18pt5_vs_Week4" "Week4_vs_2-year" "E18pt5_vs_2-year")
o="posts/post-2025-04-07.md"

{
awk '{print $1":"$2"-"$3"\tWvsE_dn\t1";}' results/2025-04-04/E18pt5_vs_Week4_diff_25p_05q_hypo.tsv 
awk '{print $1":"$2"-"$3"\tWvsE_up\t1";}' results/2025-04-04/E18pt5_vs_Week4_diff_25p_05q_hyper.tsv 
awk '{print $1":"$2"-"$3"\tYvsW_dn\t1";}' results/2025-04-04/Week4_vs_2-year_diff_25p_05q_hypo.tsv 
awk '{print $1":"$2"-"$3"\tYvsW_up\t1";}' results/2025-04-04/Week4_vs_2-year_diff_25p_05q_hyper.tsv
} | ca rc2mat - 


exit

echo '
# Differential Methylation Anlysis 
## Updates
- Increased depth by using bam files directly
- Added genomic annotations to the diff data

## Data
- fastq files: pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/250220_MY12882_fastq
- emseq results: pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/bigdata/emseq
- bedGraph files: pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/bedGraph (need to be merged by sample)

## Method
- Alignment was performed using Bismark.
- (Pair-wise) Differential analysis was conducted using methylKit.
- Bismark coverage tables were converted into methylKit input format.
- Strand-specific counts were merged, combining OriginTop (OT) and OriginBot (OB) strand frequency and counts [link]('$idir').
- Logistic regression was applied to test the odds ratio of methylation proportions between two groups using methylKit::calculateDiffMeth.

## Results
- methylKit does not support multiple comparisons.
- We can infer relationships from pairwise comparisons.

| | E18pt5 vs Week4 | Week4 vs 2-year | E18pt5 vs 2-year |
|-|-|-|-| ' > $o 


#E18pt5_vs_2-year_cluster.png                        
#E18pt5_vs_2-year_cor.png                            
#E18pt5_vs_2-year_diff.tsv                           
#E18pt5_vs_2-year_diff_25p_05q_hyper.tsv             
#E18pt5_vs_2-year_diff_25p_05q_hypo.tsv              
#E18pt5_vs_2-year_diff_25p_05q_per_chrom.png         
#E18pt5_vs_2-year_diff_25p_05q_per_cpgfeature.png    
#E18pt5_vs_2-year_diff_25p_05q_per_genomefeature.png 

printf "| Correlation |" >> $o 
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    printf "%s|" "![cor]($idir/${c}_cor.png)" >> $o
done 
printf "\n" >> $o

printf "| Cluster|" >> $o 
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    printf "%s|" "![cluster]($idir/${c}_cluster.png)" >> $o
done 
printf "\n" >> $o

printf "| Diff Sites |" >> $o 
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    i=${c}_diff_25p_05q_hyper.tsv; 
    j=${c}_diff_25p_05q_hypo.tsv; 
    printf "%s,%s|" "[hyper]($idir/$i)" "[hypo]($idir/$j)" >> $o
done 
printf "\n" >> $o

printf "| Chrom Profile|" >> $o
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    printf "!%s|" "[chrom_fig]($idir/${c}_diff_25p_05q_per_chrom.png)" >> $o
done
printf "\n" >> $o

printf "| CpG Feature |" >> $o
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    printf "!%s|" "[cpg_fig]($idir/${c}_diff_25p_05q_per_cpgfeature.png)" >> $o
done
printf "\n" >> $o

printf "| Genome Feature |" >> $o
for c in E18pt5_vs_Week4  Week4_vs_2-year  E18pt5_vs_2-year ; do
    printf "!%s|" "[genome_fig]($idir/${c}_diff_25p_05q_per_genomefeature.png)" >>$o
done
printf "\n" >> $o

exit
printf "| hypo Sites |" >> $o 
for i in ${iter[@]};do 
    ii=${i}_diff_25p_05q_hypo.tsv; printf "%s|" "[$ii]($idir/$ii)" >> $o
done 
printf "\n" >> $o

printf "| chrom Profiles |" >> $o 
for i in ${iter[@]};do 
    ii=${i}_diff_25p_05q_hypo.tsv; printf "%s|" "[$ii]($idir/$ii)" >> $o
done 
printf "\n" >> $o


exit
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

