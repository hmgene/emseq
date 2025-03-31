
input="
bigdata/methylkit/20250219_1_E18pt5.gt3.txt
bigdata/methylkit/20250219_2_E18pt5.gt3.txt
bigdata/methylkit/20250219_3_E18pt5.gt3.txt
bigdata/methylkit/20250219_4_E18pt5.gt3.txt
bigdata/methylkit/20250219_5_Week4.gt3.txt
bigdata/methylkit/20250219_6_Week4.gt3.txt
bigdata/methylkit/20250219_7_Week4.gt3.txt
bigdata/methylkit/20250219_8_Week4.gt3.txt
bigdata/methylkit/20250219_9_2-year.gt3.txt
bigdata/methylkit/20250219_10_2-year.gt3.txt
bigdata/methylkit/20250219_11_2-year.gt3.txt
bigdata/methylkit/20250219_12_2-year.gt3.txt
"

#ctr="E18pt5"
#trt="2-year"
ctr="Week4"
trt="2-year"
output=paste0("results/",ctr,"_vs_",trt)

library(methylKit)
file.list = as.list(unlist(strsplit(trimws(input), "\n")))
file.list = file.list[ grepl(ctr, file.list) | grepl(trt, file.list) ]


t=rep(0,length(file.list)); t[ grep(trt,file.list) ] = 1
s=rep("",length(file.list))
s[ grep(ctr,file.list) ] = paste0(ctr,"_",1:(sum(t==0)))
s[ grep(trt,file.list) ] = paste0(trt,"_",1:(sum(t==1)))

myobj=methRead(file.list, sample.id=as.list(s), assembly="mm10", treatment=t, context="CpG", mincov = 10 )

for( i in 1:length(myobj)){
        d=myobj[[i]]
        f=paste0("figures/cpg_freq_", d@sample.id,".png")
        png(file=f,width=800)
        getMethylationStats(d,plot=TRUE,both.strands=T)
        dev.off();

        f=paste0("figures/cpg_cov_", d@sample.id,".png")
        png(file=f,width=800)
        getCoverageStats(myobj[[2]],plot=TRUE,both.strands=T)
        dev.off()
}       

meth=unite(myobj, destrand=T,min.per.group=2L)
png(paste0(output,"_cor.png"))
getCorrelation(meth,plot=TRUE)
dev.off();
png(paste0(output,"_cluster.png"))
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off();

myDiff=calculateDiffMeth(meth)
write.table(myDiff,file=paste0(output,"_diff.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.05,type="hyper")
write.table(myDiff25p.hyper,file=paste0(output,"_diff_25p_05q_hyper.tsv"),col.names=T,row.names=F,quote=F,sep="\t")
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.05,type="hypo")
write.table(myDiff25p.hypo,file=paste0(output,"_diff_25p_05q_hypo.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

png(paste0(output,"_diff_25p_05q_per_chrom.png"))
diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.05, meth.cutoff=25)
dev.off()

library(genomation)
gene.obj=genomation::readTranscriptFeatures("data/ncbiRefseq.bed12.gz")
diffAnn.hyper=annotateWithGeneParts(as(myDiff25p.hyper,"GRanges"),gene.obj)
diffAnn.hypo=annotateWithGeneParts(as(myDiff25p.hypo,"GRanges"),gene.obj)

png(paste0(output,"_diff_25p_05q_per_genomefeature.png"))
par(mfrow=c(1,2))
plotTargetAnnotation(diffAnn.hyper,precedence=TRUE, main="hyper methylation annotation")
plotTargetAnnotation(diffAnn.hypo,precedence=TRUE, main="hypo methylation annotation")
dev.off()



