## mapping summary

input="
Y1 bigdata/methylseq/2-year_1/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y2 bigdata/methylseq/2-year_2/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y3 bigdata/methylseq/2-year_3/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y4 bigdata/methylseq/2-year_4/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E1 bigdata/methylseq/E18pt5_1/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E2 bigdata/methylseq/E18pt5_2/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E3 bigdata/methylseq/E18pt5_3/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E4 bigdata/methylseq/E18pt5_4/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W1 bigdata/methylseq/Week4_1/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W2 bigdata/methylseq/Week4_2/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W3 bigdata/methylseq/Week4_3/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W4 bigdata/methylseq/Week4_4/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
"  
echo "sample group
E1 E
E2 E
E3 E
E4 E
W1 W
W2 W
W3 W
W4 W
Y1 Y
Y2 Y
Y3 Y
Y4 Y" > meta.txt

#scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/filtered.*.*bp.anova.tsv.gz" bigdata/

#for x in 3x 5x;do
#for b in 1bp 10bp;do
#for i in E W Y;do
#    gunzip -dc bigdata/filtered.$x.$b.anova.tsv.gz | hm cutn - chrom,start,end,mean_${i} | tail -n+2 | awk -vOFS="\t" '{print $1,int($2),int($3),$4;}'  > results/mean_${i}.$x.$b.bedGraph
#done
#done
#done

#echo " 
#3x10b bigdata/filtered.3x.10bp.anova.tsv.gz
#5x10b bigdata/filtered.5x.10bp.anova.tsv.gz
#3x1b bigdata/filtered.3x.1bp.anova.tsv.gz
#5x1b bigdata/filtered.5x.1bp.anova.tsv.gz
#" | hm bismark-merge-anova -  o

#hm bismark-anno bigdata/filtered.3x.10bp.anova.tsv.gz mm10 results/filtered.3x.10bp.anova.anno
t=read.table(text="
Lin28b
Hmga2
Hic2
Igf2bp2
Igf2bp3
")
tt=fread("results/filtered.3x.10bp.anova.anno.tsv")
tt1=tt[`Gene Name` %in% t$V1] 
p= as.matrix(tt1[,grep("_perc",names(tt1),value=T),with=F] )
rownames(p) <- make.unique(paste0(tt1$`Gene Name`,".",substr(tt1$Annotation,1,10),".",tt1$chrom,":",tt1$start,"-",tt1$end))
pdf("results/target_heatmap.pdf")
Heatmap(p, row_names_gp = grid::gpar(fontsize = 5))
dev.off()
pdf("results/target_heatmap_scaled.pdf")
Heatmap(t(scale(t(p))), row_names_gp = grid::gpar(fontsize = 5))
dev.off()
###start trend
library(ggplot2)
library(data.table)
library(reshape2)

tt[, `:=` (
  trend_EW = fifelse(mean_W > mean_E, "up", fifelse(mean_W < mean_E, "dn", "no_change")),
  trend_WY = fifelse(mean_Y > mean_W, "up", fifelse(mean_Y < mean_W, "dn", "no_change"))
)]
# Now classify the full trend E → W → Y
tt[, trend_class := paste(trend_EW, trend_WY, sep = "_")]

# Melt data to long format for ggplot (E, W, Y as condition)
tt_long <- melt(tt, id.vars = c("trend_class"), measure.vars = c("mean_E", "mean_W", "mean_Y"),
                variable.name = "condition", value.name = "methylation")
# Clean up condition names
tt_long[, condition := gsub("mean_", "", condition)]
trend_counts <- tt[, .N, by = trend_class]

# Create a new label including count
trend_counts[, trend_label := paste0(trend_class, "\n(n=", N, ")")]

# Merge back with the long-form data
tt_long <- merge(tt_long, trend_counts[, .(trend_class, trend_label)], by = "trend_class")

ggplot(tt_long, aes(x = trend_label, y = methylation, fill = condition)) +
  geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.8)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2, 
               position = position_dodge(width = 0.8)) +
  labs(title = "Methylation Trends by Class", 
       x = "Trend Class (E→W→Y)", 
       y = "Methylation (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fwrite(file="results/filtered.3x.10bp.anova.anno.trend.tsv",tt)


###end




x= unlist(lapply(tt$Annotation, function(x) strsplit(x," ")[[1]][1] )
library(ggplot2)
tt[, group := sapply(Annotation, function(x) strsplit(x, " ")[[1]][1])]
mean_vals <- tt[, .(
  mean_E = mean(mean_E, na.rm = TRUE),
  mean_W = mean(mean_W, na.rm = TRUE),
  mean_Y = mean(mean_Y, na.rm = TRUE)
), by = group]

mean_vals_long <- melt(mean_vals, id.vars = "group", variable.name = "condition", value.name = "mean_value")
ggplot(mean_vals_long, aes(x = group, y = mean_value, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + labs(title = "Mean Values per Group", x = "Group", y = "Mean Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## distributions
tt_melt <- melt(tt, id.vars = "group",
                measure.vars = c("mean_E", "mean_W", "mean_Y"),
                variable.name = "condition",
                value.name = "methylation")

# Remove NA values
tt_melt <- tt_melt[!is.na(methylation)]

# Plot: Methylation distribution for each group across 3 conditions
ggplot(tt_melt, aes(x = methylation, fill = condition)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~group, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of Methylation by Group (E, W, Y)",
       x = "Methylation (%)",
       y = "Density") +
  theme(text = element_text(size = 12))




## run this on HPC
##echo "$input"| grep -v "^$" | hm bismark-merge-cov - 5 | gzip -c > bigdata/merged.x5.cov.gz 
##scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/merged.x5.cov.gz bigdata/
##scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/data/emseq/merged.3x.cov.gz bigdata/
##hm bismark-anova bigdata/filtered.3x.cov.gz meta.txt 0.05 > bigdata/filtered.cov.anova.tsv
exit

#hm bismark-cov-bin bigdata/filtered.cov.gz 1 | hm bismark-anova -  meta.txt 0.001 > bigdata/filtered.cov.${bz}bp.anova.p-3.tsv
#hm bismark-cov-groupfilter bigdata/merged.5x.cov.gz meta.txt 2 5 | gzip -c > bigdata/filtered.5x.cov.gz 
#hm bismark-cov-groupfilter bigdata/merged.5x.cov.gz meta.txt 2 5 | gzip -c > bigdata/filtered.3x.cov.gz 

exit
#hm bismark-cov-bin bigdata/filtered.cov.gz 10 | hm bismark-anova - meta.txt 0.05 > bigdata/filtered.cov.10b.anova.tsv
#for i in E W Y;do
#    hm cutn bigdata/filtered.cov.10b.anova.tsv chrom,start,end,mean_${i} | tail -n+2 | awk -vOFS="\t" '{print $1,int($2),int($3),$4;}'  > results/mean_${i}.bedGraph
#done

exit
summary(){
fir=1
for f in bigdata/methylseq/*/out/bismark/summary/bismark_summary_report.txt;do
    s=`cut -d"/" -f 3  $f`;
    if [ "$fir" -eq 0 ];then head -n 1 $f > summary.tsv; fir=0; fi
    tail -n+2 $f   >> summary.tsv
done
}

targetting(){
target=( 
Lin28b
Hmga2
Hic2
Igf2bp2
Igf2bp3
)

hm ucsc-refflat mm10 | grep -wf <( echo ${target[@]} | tr " " "\n" ) | hm ucsc-refflat2bed12 - | cut -f 1-6  | tee target.bed |\
tee >( hm bedw - 100000 > target_w100k.bed ) | hm bed5p - | hm bedw - 4000 > target_tss4k.bed

fn(){
    s=`echo $1 | cut -d"/" -f 3`
    b=10
    o=bigdata/methylseq/bw/$s.${b}bp.bw
    mkdir -p  ${o%/*};
    echo "$s.."
    gunzip -dc $1 | awk '$5+$6 > 2' | hm bismark-cov2bg - $b | hm bg2bw - mm10 > $o
};export -f fn;
parallel -j 3 fn {}  ::: bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz 
}

for f in bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz;do
    s=`echo $f | cut -d"/" -f 3`
    echo "$s $f"
done #|  bismark-merge-cov - | head

#chr start   end strand  coverage    numCs   numTs
#chr1    3053156 3053156 -   10  9   1
exit

hm mycat  bigdata/methylseq/2-year_1/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz 
exit

f="bigdata/bams/20250219_10_2-year_2_MY12882_S10_L006.deduplicated.sorted.bam" 
my.methRaw=processBismarkAln( location = f, sample.id="test1", assembly="mm10", read.context="CpG", save.folder=getwd())




. src/bed.sh; ref2bed12 data/ncbiRefSeq.txt  >  data/ncbiRefSeq.bed12



nene(){
    cat data/ncbiRefSeq.txt |\
    perl -ne 'chomp;my@d=split/\t/,$_; print join("\t",$d[2],$d[4],$d[5],$d[1],$d[12],$d[3]),"\n";' |\
    awk '!($1~/_/) && ($4~/NM_/) { $4=$5;$5=0; print $0;}' | sort -u >  data/gene.bed  
}
bw(){
    o=${1%.bedGraph.gz}.bw
    local tmpd=`mktemp -d`
    gunzip -dc $1 | tail -n+2 | sort -k1,1 -k2,3n > $tmpd/a    
    `ca home`/bin/bedGraphToBigWig $tmpd/a data/mm10.chrom.sizes $o
};export -f bw;

#gene
#parallel bw {} ::: data/20250219_*.gz
i=(
data/20250219_10_2-year_2_MY12882_S10_L006.bw	data/20250219_12_2-year_4_MY12882_S12_L007.bw	data/20250219_3_E18pt5_3_MY12882_S3_L006.bw
data/20250219_10_2-year_2_MY12882_S10_L007.bw	data/20250219_1_E18pt5_1_MY12882_S1_L006.bw	data/20250219_3_E18pt5_3_MY12882_S3_L007.bw
data/20250219_11_2-year_3_MY12882_S11_L006.bw	data/20250219_1_E18pt5_1_MY12882_S1_L007.bw	data/20250219_4_E18pt5_4_MY12882_S4_L006.bw
data/20250219_11_2-year_3_MY12882_S11_L007.bw	data/20250219_2_E18pt5_2_MY12882_S2_L006.bw	data/20250219_4_E18pt5_4_MY12882_S4_L007.bw
data/20250219_12_2-year_4_MY12882_S12_L006.bw	data/20250219_2_E18pt5_2_MY12882_S2_L007.bw
)
#computeMatrix reference-point -S ${i[@]} -R data/gene.bed  -a 3000 -b 3000 -o bigdata/gene.gz

mk-bmcov(){

bismark2methylkit(){
    echo "chrBase chr base    strand  coverage    freqC   freqT" | sed "s/  */\t/g"
    {
        #LH00547:89:22W5LCLT3:6:1101:19335:1080_1:N:0:CTCGAAAT+CTTNCCTG	+	chr5	75583394	Z
        gunzip -dc $1 | awk '{print $3,$4,"F",$5;}'
        gunzip -dc ${1/OB/} | awk '{print $3,$4,"R",$5}'
    } | perl -e 'use strict; 
        sub get{ my ($x)=@_; return defined $x ? $x : 0; }
        my %r=();
        while (<>) { chomp;  my ($c, $s, $t, $z) = split /\s+/;  
            $r{"$c\t$s\t$t"}{$z}++;  
        }
        map {my$k=$_;my($c,$s,$t)=split/\t/,$_; 
            my ($x,$y) = (get($r{$k}{"Z"}), get($r{$k}{"z"}) );
            #chrBase chr base    strand  coverage    freqC   freqT
            #chr21.9764539   chr21   9764539 R   12  25.00   75.00
            if($x+$y > 0){
                print join("\t","$c.$s",$c,$s,$t,$x+$y,
                map{ sprintf("%.2f",$_*100) } ( $x/($x+$y),$y/($x+$y))),"\n";  
            }
        } keys %r;
    '
}
    s=${1#*/bigdata/};s=${s%/bismark*}
    o=bigdata/methylkit/$s.txt
    mkdir -p ${o%/*}
    bismark2methylkit $1 > $o
};export -f mk-bmcov

#parallel mk-bmcov {} ::: \
#/Volumes/T7/git/emseq/bigdata/*/bismark/methylation_calls/methylation_calls/CpG_OT_*.deduplicated.txt.gz

## methylkit

library(methylKit)

file.list<-list.files(path = "bigdata/methylkit", pattern = "*\\.txt$", full.names = TRUE)
sample.id = as.list(basename(file.list)
treatment=rep(0,length(sample.id))
treatment[1:10]=1
myobj=methRead(as.list(file.list),sample.id=sample.id, assembly="mm10", treatment=treatment, context="CpG", mincov = 10)

meth=unite(myobj, destrand=FALSE)


