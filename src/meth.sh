bisbam2meth2(){
    #ref: methylKit/src/methCall.cpp
echo "
1_R1/1  67  5 103172224 255 40M = 103172417 233 AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:4  XX:Z:4T1T24TT7  XM:Z:....h.h........................hh....... XR:Z:CT XG:Z:CT
1_R1/2  131 5 103172417 255 40M = 103172224 -233  TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:6  XX:Z:T5T1T9T9T7T3 XM:Z:h.....h.h.........h.........h.......h... XR:Z:GA XG:Z:CT
HWI-ST986_0098:1:1101:18264:11272#0/1 0 chr1  497 255 50M * 0 0 TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA  CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI  NM:i:13 XX:Z:C4C3CC9C4C1C2CC2C6CC1C5  XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ.... XR:Z:CT XG:Z:CT
" | grep -v "^$" | sed -E "s/ +/\t/g"
}

#bisbam2meth2

bisbam2meth(){
usage="$FUNCNAME <bismark.bam> <name> <out.tsv>";
if [ $# -lt 3 ];then echo "$usage";return;fi 

    R --no-save -e '
    library(methylKit)
    o=processBismarkAln( location = "'$1'", sample.id="'$2'", assembly="mm10", read.context="CpG") #, save.folder=getwd())
    write.table(file="'$3'",o,col.names=T,row.names=F,quote=F,sep="\t")
    '
}
bisbam2meth-test(){
    bis2meth data/sample.bam sample data/sample.tsv
}
