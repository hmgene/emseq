make-bgs(){
Rscript -e '
    ## read RNA and make bedGraph for RNA-seq
    library(data.table)
    tt=fread("Beaudin.txt",skip=2) ## fl vs adult 
    tt=tt[padj<0.05,.(Name,log2FoldChange)]
    setnames(tt,old=c("log2FoldChange"),new=c("l2_FLvsAD"));
    tt[, l2_FLvsAD := mean(l2_FLvsAD,na.rm=T), by=Name]

    tmp=fread("Tan.txt",skip=2) ## adult vs fl
    tmp=tmp[padj<0.05,.(Name,log2FoldChange)]
    setnames(tmp,old=c("log2FoldChange"),new=c("l2_ADvsFL"));
    tmp[, l2_ADvsFL := mean(l2_ADvsFL,na.rm=T), by=Name]

    tt=merge(tt,tmp,all=T)
    r=fread("~/git/hmtools/data/ucsc/mm10/refFlat.txt.gz")
    r=r[,c(1,3,5,6)]
    setnames(r,c("Name","chrom","start","end"))
    tt=merge(tt,r)
    tt = tt[, .( chrom = unique(chrom), start = min(start), end = max(end), l2_FLvsAD = mean(l2_FLvsAD), l2_ADvsFL =mean(l2_ADvsFL)), by = Name]

    ## RNA to bedGraph
    res <- tt[, {
      # 1. Create IRanges object for this chromosome
      ir <- IRanges(start, end)
      disj <- disjoin(ir)   # split into non-overlapping intervals
      
      # 2. Find overlaps between disjoint intervals and original ranges
      ol <- findOverlaps(disj, ir)
      
      # 3. Split values by disjoint intervals
      l2_FL_list <- splitAsList(l2_FLvsAD[subjectHits(ol)], queryHits(ol))
      l2_AD_list <- splitAsList(l2_ADvsFL[subjectHits(ol)], queryHits(ol))
      name_list  <- splitAsList(Name[subjectHits(ol)], queryHits(ol))
      
      # 4. Build final data.table
      data.table(
        start = start(disj),
        end   = end(disj),
        l2_FLvsAD = sapply(l2_FL_list, mean, na.rm=TRUE),
        l2_ADvsFL = sapply(l2_AD_list, mean, na.rm=TRUE),
        Name = sapply(name_list, function(x) paste(unique(x), collapse=";"))
      )
      
    }, by = chrom]

    fwrite(res[,.(chrom,start,end,l2_FLvsAD)],"RNA_l2_FLvsAD_Beaudin.bedGraph",sep="\t",col.names=F)
    fwrite(res[,.(chrom,start,end,l2_ADvsFL)],"RNA_l2_ADvsFL_Tan.bedGraph",sep="\t",col.names=F)

    ## read ATAC and make bedGraph
    input=c("BM_HSC_on_FL_HSC_ATAC.xls", "FL_HSC_on_BM_HSC_ATAC.xls")
    output=c("ATAC_BMvsFL.bedGraph","ATAC_FLvsBM.bedGraph");
    for( i in 1:length(input)){
        tt=fread(input[i])
        fwrite(tt[,.(chr,start,end,fold_enrichment)],output[i],sep="\t",col.names=F)
        
    }

    tt=fread("BM_HSC_on_FL_HSC_ATAC.xls")[,.(chr,start,end,fold_enrichment)]
    setnames(tt,old="fold_enrichment","ATAC_BMvsFL")
    fwrite(tt[,.(chr,start,end,fold_enrichment)],
	tmp=fread("FL_HSC_on_BM_HSC_ATAC.xls")[,.(chr,start,end,fold_enrichment)]
    setnames(tmp,old="fold_enrichment","ATAC_FLvsBM")
    tt=merge(tt,tmp)

'
}


Rscript -e 'library(rmarkdown);render("README.Rmd")'
git add -A
git commit -am "rmd"
git push
