

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


