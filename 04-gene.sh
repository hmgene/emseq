
tail -n+2 results/summary_1k_diffmeth.tsv |\
cut -f 1 | sed "s/[.-]/\t/g" | windowBed -a stdin -b data/gene.bed -w 10000 |\
awk -v OFS="\t" '{ print $1"."$2"-"$3, $7;}'  | sort -u > results/summary_1k_diffmeth.tsv.gene



R --no-save -e '

tt=read.table("results/summary_1k_diffmeth.tsv",sep="\t",header=T)
x=read.table("results/summary_1k_diffmeth.tsv.gene",sep="\t",header=F)
colnames(x)=c("rid","gene")
d = merge(tt, x, by = "rid", all.x = TRUE)
d$rid = paste0(d$rid,ifelse( is.na(d$gene),"",paste0(":",d$gene)))

wrtie.table(d,file="results/summary_1k_diffmeth_gene.tsv",sep="\t",col.names=T,row.names=F,quote=F,sep="\t)



