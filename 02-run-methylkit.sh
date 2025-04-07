

input="
bigdata/meth/10_2-year_2.tsv  bigdata/meth/12_2-year_4.tsv  bigdata/meth/2_E18pt5_2.tsv  bigdata/meth/4_E18pt5_4.tsv  bigdata/meth/6_Week4_2.tsv  bigdata/meth/8_Week4_4.tsv
bigdata/meth/11_2-year_3.tsv  bigdata/meth/1_E18pt5_1.tsv   bigdata/meth/3_E18pt5_3.tsv  bigdata/meth/5_Week4_1.tsv   bigdata/meth/7_Week4_3.tsv  bigdata/meth/9_2-year_1.tsv
"
odir="post/2025-04-04"
ctr="E18pt5";trt="Week4"
#ctr="E18pt5";trt="2-year"
#ctr="Week4";trt="2-year"
output=paste0(odir,"/",ctr,"_vs_",trt)


library(methylKit)
library(genomation)

file.list = as.list(unlist(strsplit(trimws(input), "\\s+",perl=T)))
file.list = file.list[ grepl(ctr, file.list) | grepl(trt, file.list) ]

t=rep(0,length(file.list)); t[ grep(trt,file.list) ] = 1
s=rep("",length(file.list))
s[ grep(ctr,file.list) ] = paste0(ctr,"_",1:(sum(t==0)))
s[ grep(trt,file.list) ] = paste0(trt,"_",1:(sum(t==1)))

myobj <- mapply(function(f, i) {
    df=read.table(f,header=T)
    new("methylRaw", df, sample.id = i, assembly = "mm10", context = "CpG", resolution = "base")
}, file.list , s, SIMPLIFY = FALSE)
myobj=methylRawList(myobj,treatment=t)

cpg.obj=readFeatureFlank( "data/cpgIslandExt.bed", feature.flank.name=c("CpGi","shores"))
gene.obj=genomation::readTranscriptFeatures("data/ncbiRefseq.bed12.gz")


for( i in 1:length(myobj)){
        d=myobj[[i]]
        f=paste0(gsub("/[^/]+$","",output,perl=T),"/",d@sample.id,"_cpg_freq.png")

        png(file=f,width=800)
        getMethylationStats(d,plot=TRUE,both.strands=T)
        dev.off();

        f=paste0(gsub("/[^/]+$","",output,perl=T),"/",d@sample.id,"_cov_freq.png")
        png(file=f,width=800)
        getCoverageStats(d,plot=TRUE,both.strands=T)
        dev.off()
}       

meth=unite(myobj, destrand=T,min.per.group=2L)

png(paste0(output,"_cor.png"))
getCorrelation(meth,plot=TRUE)
dev.off();

png(paste0(output,"_cluster.png"))
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off();


## differential 
myDiff=calculateDiffMeth(meth)
x=as.data.frame(myDiff@.Data,col.names=myDiff@names)
y=annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)
z=annotateWithFeatureFlank(as(myDiff,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
anno=cbind(as.data.frame(myDiff@.Data,col.names=myDiff@names), y@dist.to.TSS,y@members,z@members)
write.table(anno,file=paste0(output,"_diff.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.05,type="hyper")
#write.table(myDiff25p.hyper,file=paste0(output,"_diff_25p_05q_hyper.tsv"),col.names=T,row.names=F,quote=F,sep="\t")
write.table(anno[anno$qvalue <= 0.05 & anno$meth.diff >= 0.25,],file=paste0(output,"_diff_25p_05q_hyper.tsv"),col.names=T,row.names=F,quote=F,sep="\t")


myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.05,type="hypo")
#write.table(myDiff25p.hypo,file=paste0(output,"_diff_25p_05q_hypo.tsv"),col.names=T,row.names=F,quote=F,sep="\t")
write.table(anno[anno$qvalue <= 0.05 & anno$meth.diff <= -0.25,],file=paste0(output,"_diff_25p_05q_hypo.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

png(paste0(output,"_diff_25p_05q_per_chrom.png"))
diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.05, meth.cutoff=25)
dev.off()

diffAnn.hyper=annotateWithGeneParts(as(myDiff25p.hyper,"GRanges"),gene.obj)
diffAnn.hypo=annotateWithGeneParts(as(myDiff25p.hypo,"GRanges"),gene.obj)

png(paste0(output,"_diff_25p_05q_per_genomefeature.png"))
par(mfrow=c(1,2))
plotTargetAnnotation(diffAnn.hyper,precedence=TRUE, main="hyper methylation annotation")
plotTargetAnnotation(diffAnn.hypo,precedence=TRUE, main="hypo methylation annotation")
dev.off()


cpg.obj=readFeatureFlank( "data/cpgIslandExt.bed", feature.flank.name=c("CpGi","shores"))
diffCpGann.hyper=annotateWithFeatureFlank(as(myDiff25p.hyper,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
diffCpGann.hypo=annotateWithFeatureFlank(as(myDiff25p.hypo,"GRanges"), cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")

png(paste0(output,"_diff_25p_05q_per_cpgfeature.png"))
par(mfrow=c(1,2))
plotTargetAnnotation(diffCpGann.hypo,col=c("green","gray","white"), main="hypo methylation annotation")
plotTargetAnnotation(diffCpGann.hyper,col=c("green","gray","white"), main="hyper methylation annotation")
dev.off();



