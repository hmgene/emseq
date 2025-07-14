odir=results/2025-07-14
mkdir -p $odir


echo '
library(data.table)
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
tt1 = tt[tt[, apply(.SD < 0.00001, 1, any), .SDcols = pv_cols]]
#pt=0.001; tt1 = tt[ rna_fl_pv < pt & atac_fl_pv < pt & rna_bm_pv < pt & atac_bm_pv < pt ]

j=setdiff(names(tt1),c("gene",pv_cols));1
m=as.matrix(tt1[,..j]) #(atac_bm,rna_bm)]
row.names(m) = tt1$gene

#d=fread("bigdata/filtered.3x.10bp.anova.tsv.gz")
d=fread("results/2025-06-03-3x-depth/filtered.3x.10bp.anova.anno.trend.tsv.gz")
d1=d[,c("Gene Name",grep("_perc",names(d),value=T)),with=F]
setnames(d1,"Gene Name","gene")
d1 = d1[, lapply(.SD, mean, na.rm = TRUE), by = gene]
d1=d1[gene %in% tt1$gene,]


m1=as.matrix(d1[, grep("E|W",names(d1),value=T),with=F])
m1[is.na(m1)]=0
row.names(m1)=d1$gene
m1=m1[rownames(m),]


M=scale(m)
M1=t(scale(t(m1)))

row_order <- hclust(dist(M))  # or dist(m1) if appropriate
row_dend <- as.dendrogram(row_order)

h1 <- Heatmap(M, name = "Scaled Signal", cluster_rows = row_dend, show_row_names = FALSE, column_title = "ATAC / RNA")
h2 <- Heatmap(M1, name = "Percent", cluster_rows = row_dend, show_row_names = FALSE, column_title = "Methylation")
draw(h2 + h1)

Heatmap(cbind(M,M1))




data/from_pub/FL_HSC_on_BM_HSC_ATAC.xls.anno
atac_bm_fl=fread("results/2025-07-14/bm_vs_fl.anno")
atac_fl_bm=fread("results/2025-07-14/fl_vs_bm.anno")

tt = rna_fl_bm[, .(gene = Name, rna_fl = log2FoldChange)]
tt = merge(tt,rna_bm_fl[, .(gene = Name, rna_bm = log2FoldChange)])1
tt = merge(tt,rna_bm_fl[, .(gene = Name, rna_bm = log2FoldChange)])1





'



