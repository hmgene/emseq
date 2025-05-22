## mapping summary

echo "File Total_Reads Aligned_Reads" | tr " " "\t" > summary.tsv
for f in bigdata/methylseq/*/out/bismark/summary/bismark_summary_report.txt;do
    s=`cut -d"/" -f 3  $f`;
    tail -n+2 $f |  cut -f 1-3 >> summary.tsv
done

exit

target=( 
Lin28b
Hmga2
Hic2
Igf2bp2
Igf2bp3
)

#hm ucsc-refflat mm10 | grep -wf <( echo ${target[@]} | tr " " "\n" ) | hm ucsc-refflat2bed12 - | cut -f 1-6  | hm bed5p - | hm bedw - 4000 > target_tss4k.bed

fn(){
    s=`echo $1 | cut -d"/" -f 3`
    b=10
    o=bigdata/methylseq/bw/$s.${b}bp.bw
    mkdir -p  ${o%/*};
    echo "$s.."
    gunzip -dc $1 | awk '$5+$6 > 2' | hm bismark-cov2bg - $b | hm bg2bw - mm10 > $o
};export -f fn;
#parallel -j 3 fn {}  ::: bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz 

for f in bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*.deduplicated.bismark.cov.gz;do
    s=`echo $f | cut -d"/" -f 3`
    echo "$s $f"
done > input.txt 

import pandas as pd
import sys
inp = pd.read_csv(sys.stdin,sep=" ",header=None)
d=null
for i, row in inp.iterrows():
    s = row[0]
    f = row[1]
    tt=pd.read_csv(f,header=None,sep="\t")
    tt.columns=["chrom","start","end"]+[f"{i}_{s}"  for i in ["perc","numC","numT"]] 
    tt = tt.loc[tt[f"numC_{s}"] + tt[f"numT_{s}"] > 2]
    if d is None:
        d=tt
    else:
        d=d.merge(tt,on=["chrom", "start", "end"], how="outer")

d.to_csv(sys.stdout,sep="\t",index=False)

 Rscript -e 'tt=read.table("input.txt",header=F)
    library(data.table)
    d=NULL;
    for( i in 1:nrow(tt)){
        if(is.null(d)){
            d=fread(tt[i,2], col.names=c("chrom","start","end",paste(tt[i,1],c("perc","numC","numT"),sep="_")))
        }else{
            d=merge(d,fread(tt[i,2], col.names=c("chrom","start","end",paste(tt[i,1],c("perc","numC","numT"),sep="_"))),all=T)
        }
    }
'

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


