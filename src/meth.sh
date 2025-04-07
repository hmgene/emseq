
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
