odir="results/2025-07-14"
mkdir -p $odir


echo '
library(data.table)
library(ComplexHeatmap)
#FL versus adult HSCs
tt=fread("data/from_pub/Beaudin.txt",skip=2)
tt=tt[,.(gene=Name, rna_fl = log2FoldChange, rna_fl_pv = pvalue)]
tt=merge(tt,fread("data/from_pub/Tan.txt",skip=2)[,.(gene=Name, rna_bm = log2FoldChange, rna_bm_pv=pvalue)])

a <- fread("data/from_pub/BM_HSC_on_FL_HSC_ATAC.xls")[ , .(Symbol = unlist(strsplit(Symbol, ";", fixed = TRUE)), fold_enrichment, `-log10(pvalue)`), by = .(chr, start, end) ]
a <- a[!is.na(Symbol) & Symbol != "", .( atac_bm = max(fold_enrichment, na.rm = TRUE), atac_bm_pv = min(10^(-`-log10(pvalue)`), na.rm = TRUE)), by = .(gene = Symbol)]
b <- fread("data/from_pub/FL_HSC_on_BM_HSC_ATAC.xls")[ , .(Symbol = unlist(strsplit(Symbol, ";", fixed = TRUE)), fold_enrichment, `-log10(pvalue)`), by = .(chr, start, end) ]
b <- b[!is.na(Symbol) & Symbol != "", .( atac_fl = max(fold_enrichment, na.rm = TRUE), atac_fl_pv = min(10^(-`-log10(pvalue)`), na.rm = TRUE)), by = .(gene = Symbol)]
tt=merge(tt,a)
tt=merge(tt,b)

pv_cols <- grep("_pv$", names(tt), value = TRUE)
#tt1 = tt[tt[, apply(.SD < 0.00001, 1, any), .SDcols = pv_cols]]
pt=0.001; tt1 = tt[ rna_fl_pv < pt & atac_fl_pv < pt & rna_bm_pv < pt & atac_bm_pv < pt ]

j=setdiff(names(tt1),c("gene",pv_cols));1
m=as.matrix(tt1[,..j]) #(atac_bm,rna_bm)]
row.names(m) = tt1$gene

#d=fread("bigdata/filtered.3x.10bp.anova.tsv.gz")
d=fread("results/2025-06-03-3x-depth/filtered.3x.10bp.anova.anno.trend.tsv.gz")
for( j in c("mean_E","mean_W","mean_Y")){
    fwrite(file=paste0(odir,"/",j,".bedGraph.gz"),sep="\t",col.names=F, d[,.(chrom,start,end,value=get(j))])
}
fwrite(file=paste0(odir,"/EtoW_up.bed"), d[ trend_EW == "up",.(chrom,start,end)],col.names=F,sep="\t")
fwrite(file=paste0(odir,"/EtoW_dn.bed"), d[ trend_EW == "dn",.(chrom,start,end)],col.names=F,sep="\t")
fwrite(file=paste0(odir,"/WtoY_up.bed"), d[ trend_WY == "up",.(chrom,start,end)],col.names=F,sep="\t")
fwrite(file=paste0(odir,"/WtoY_dn.bed"), d[ trend_WY == "dn",.(chrom,start,end)],col.names=F,sep="\t")



d1=d[,c("Gene Name",grep("_perc",names(d),value=T)),with=F]
setnames(d1,"Gene Name","gene")
d1 = d1[, lapply(.SD, mean, na.rm = TRUE), by = gene]
d1=d1[gene %in% tt1$gene,]
m1=as.matrix(d1[, grep("E|W",names(d1),value=T),with=F])
m1[is.na(m1)]=0
row.names(m1)=d1$gene

common_genes <- intersect(rownames(m), rownames(m1))
m_common <- m[common_genes, ]
m1_common <- m1[common_genes, ]
M <- scale(m_common)
M1 <- t(scale(t(m1_common)))
row_order <- hclust(dist(M))
row_dend <- as.dendrogram(row_order)

library(ComplexHeatmap)
h1 <- Heatmap(M, name = "Scaled Signal", cluster_rows = row_dend,
              show_row_names = FALSE, column_title = "ATAC / RNA")
h2 <- Heatmap(M1, name = "Percent", cluster_rows = row_dend,
              show_row_names = FALSE, column_title = "Methylation")
odir="results/2025-07-14";
mm = cbind(M, M1); mm = as.data.table(mm, keep.rownames = "gene") 

fwrite(file=paste0(odir,"/meth_rna_atac.csv"), mm)

pdf(file=paste0(odir,"/meth_vs_rnaatac.heatmap.pdf"))
draw(h2 + h1)
dev.off()

pdf(file=paste0(odir,"/rnaatac_vs_meth.heatmap.pdf"))
draw(h1 + h2)
dev.off()

h=Heatmap(cbind(M,M1))
pdf(file=paste0(odir,"/meth_rna_atac_mix.heatmap.pdf"))
draw(h)
dev.off();


h1 <- Heatmap(M, name = "Scaled Signal", cluster_rows = row_dend,

#data/from_pub/FL_HSC_on_BM_HSC_ATAC.xls.anno
#atac_bm_fl=fread("results/2025-07-14/bm_vs_fl.anno")
#atac_fl_bm=fread("results/2025-07-14/fl_vs_bm.anno")
#
#tt = rna_fl_bm[, .(gene = Name, rna_fl = log2FoldChange)]
#tt = merge(tt,rna_bm_fl[, .(gene = Name, rna_bm = log2FoldChange)])1
#tt = merge(tt,rna_bm_fl[, .(gene = Name, rna_bm = log2FoldChange)])1





'



