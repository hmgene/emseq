input="bigdata/2025-04-11/merged.tsv"
output="results/2025-05-08/anova"
odir=${output%/*}
mkdir -p $odir 
o=$odir/README.md
. src/util.sh;
cutn ${output}_multi.tsv chr,start,end,min_pval,gW_10,gY_10 | awk -v OFS="\t" '$4<1{ b=int($2/10)*10;print $1,b,b+10,$5;}' | sort -k1,1 -k2,3n > gW.bedGraph 
cutn ${output}_multi.tsv chr,start,end,min_pval,gW_10,gY_10 | awk -v OFS="\t" '$4<1{ b=int($2/10)*10;print $1,b,b+10,$6;}' | sort -k1,1 -k2,3n > gY.bedGraph 
exit

mamba activate r4.4
Rscript <( echo '
#https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
    input="results/2025-05-08/anova_multi.tsv"
    output="results/2025-05-08/anova_multi"
    library(ChIPseeker)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene      
    library(clusterProfiler)
    library(data.table)
    d=fread("results/2025-05-08/anova_multi.tsv")
    p=d[min_pval<0.05,c("chr","start","end","gW_10","gY_10")]
    a=d[,c("chr","start","end","gW_10","gY_10")]

d_gr = GRanges( seqnames = d$chr, ranges = IRanges(d$start,d$end), strand=d$strand)
b <- fread("target_region.bed", col.names = c("chr", "start", "end","gene","score","strand"))
b_gr<- GRanges(seqnames = b$chr, ranges = IRanges(start = b$start + 1, end = b$end))  # BED is 0-based

g=d[queryHits(findOverlaps(d_gr,b_gr)),]


    p1 = GRanges( seqnames = p$chr, ranges = IRanges(p$start,p$end), strand=p$strand,score=p$gW_10)
    p2 = GRanges( seqnames = p$chr, ranges = IRanges(p$start,p$end), strand=p$strand,score=p$gY_10)
    p12 <- list(gW_10 = p1, gY_10 = p2)
    a1 = GRanges( seqnames = a$chr, ranges = IRanges(a$start,a$end), strand=a$strand,score=a$gW_10)
    a2 = GRanges( seqnames = a$chr, ranges = IRanges(a$start,a$end), strand=a$strand,score=a$gY_10)
    a12 = list(gW_10 = a1, gY_10 = a2)
    
    g1 = GRanges( seqnames = g$chr, ranges = IRanges(g$start,g$end), strand=g$strand,score=g$gW_10)
    g2 = GRanges( seqnames = g$chr, ranges = IRanges(g$start,g$end), strand=g$strand,score=g$gY_10)
    g12 = list(gW_10 = g1, gY_10 = g2)
plotPeakProf2(g12, upstream = 10000, downstream = 10000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)

plotPeakProf2(g12, upstream = 100000, downstream = 100000,
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

plotPeakProf2(a12, upstream = 10000, downstream = 10000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)


    png(paste0(output,"_chromview.png"))
    covplot(p12, weightCol="score",title="sig CpGMeth Week and Year / E",fill_color = c("red","blue"))
    dev.off()



# Install if needed
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)

# Example: list of Entrez gene IDs or gene symbols
genes <- c("Hmga2", "Igf2bp2", "Igf2bp3", "Lin28b", "Hic2")

# If using gene SYMBOLs, convert to Entrez IDs:
gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Run GO enrichment
go_result <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",           # "BP"=Biological Process, "MF", "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# View results
head(go_result)


peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000)

txdb1 <- transcripts(txdb)

target_genes=( Hic2 Lin28b Igf2bp2 Igf2bp3 and Hmga2 )


region_list <- list(geneX = txdb1, geneY = txdb2)
peakHeatmap_multiple_Sets(peak = files[[4]],
                          upstream = 1000,downstream = 1000,
                          by = c("geneX","geneY"),
                          type = "start_site",
                          TxDb = region_list,nbin = 800)


')




target_genes=( Hic2 Lin28b Igf2bp2 Igf2bp3 and Hmga2 )
grep -wf <( echo ${target_genes[@]} | tr " " "\n" ) data/gene.bed | sort -k1,1 -k2,3n |\
tee >( awk -v OFS="\t" '{ print $1, $2-10000000,$3+10000000,$4,$5,$6;}' > target_region.bed) |\
closestBed -d -a stdin -b <( cutn ${output}_multi.tsv chr,start,end,min_pval | sort -k1,1 -k2,3n )

cat data/gene.bed | sort -k1,1 -k2,3n |\
windowBed -w 10000 -a stdin -b <( cutn ${output}_multi.tsv chr,start,end,min_pval | awk '$4<0.05' | sort -k1,1 -k2,3n )

exit


. ./src/stat.sh;

run-anova(){
    for b in 1 10 50;do
        anova $input $b 1 > ${output}_${b}bp.tsv
    done
}
run-mult(){
    multibin ${output}_multi results/2025-05-08/anova_10bp.tsv	results/2025-05-08/anova_1bp.tsv	results/2025-05-08/anova_50bp.tsv
}



exit



library(cluster)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)



library(dplyr)
library(tidyr)
library(stringr)
library(cluster)

# Step 1: Convert the matrix to a long format with sample metadata
df_long <- as.data.frame(m) %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "value") %>%
  mutate(
    group = str_extract(sample, "^[^_]+"),                   # e.g., "E18pt5"
    replicate = str_extract(sample, "(?<=_)[0-9]+(?=_)"),    # e.g., "1"
    bin_size = str_extract(sample, "perc_\\d+")              # e.g., "perc_50"
  )


results <- df_wide %>%
  group_by(group, bin_size) %>% summarise()
  summarize(silhouette_score = compute_sil(cur_data()), .groups = "drop")

print(results)




# Example: convert your matrix to a data frame (if needed)
# df <- as.data.frame(your_matrix)

# Get long-form data
df_long <- as.data.frame(m) %>%
  rownames_to_column("feature_id") %>%  # if rownames are meaningful
  pivot_longer(-feature_id, names_to = "sample", values_to = "value") %>%
  mutate(
    group = str_extract(sample, "^[^_]+"),
    replicate = str_extract(sample, "(?<=_)[0-9]+(?=_)"),
    bin_size = str_extract(sample, "perc_\\d+")
  )

# Now pivot wider so each sample is a row, each feature is a column
df_wide <- df_long %>%
  pivot_wider(names_from = feature_id, values_from = value)

# Compute silhouette score for each group+bin+replicate combo
compute_sil <- function(df_group) {
  df_feat <- df_group %>% select(-sample, -group, -replicate, -bin_size)
  if (nrow(df_feat) < 2) return(NA)  # need at least 2 samples
  # Dummy clusters: assume each replicate is its own cluster
  # OR: apply real clustering (e.g., kmeans) if you have ground truth
  km <- kmeans(df_feat, centers = 2)
  d <- dist(df_feat)
  s <- silhouette(km$cluster, d)
  mean(s[, "sil_width"])
}

results <- df_wide %>%
  group_by(group, bin_size) %>%
  summarise(silhouette_score = compute_sil(cur_data()), .groups = "drop")

print(results)

    
')

       

# Clean and rename
result <- tt_joined[, .(chr, start = i.start, end = i.end, pval_1 = pval_10, pval_10 = i.pval_10)]
    
    tt = foverlaps(tt,tt1, by.x = c("chr", "start", "end"), type = "within", nomatch = NA)


')

exit




    pvals[is.na(pvals)]=1

    i= pvals < 0.05
    df=cbind(tt[,1:4],p,pvals)
    write.table(df,file=paste0(o,"_anova.tsv"),sep="\t",row.names=F,col.names=T,quote=F)

    library(ggplot2)
    library(reshape2)
    p_clean <- t(scale(t(p[i,])))  # or use: na.omit(p)
    p_clean$id <- 1:nrow(p_clean)
    meth_long <- melt(p_clean, id.vars = "id", variable.name = "Sample", value.name = "Methylation")
    ggplot(meth_long, aes(x = Sample, y = Methylation, fill = Sample)) +
      geom_violin(trim = FALSE, scale = "width") +
      theme_minimal() +
      labs(title = "Methylation Distribution per Sample", x = "Sample", y = "Methylation Level") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ## clustyering
    m=p[i,]
    m=t(scale(t(m)))
    h =Heatmap(m)
    row_dend <- row_dend(h)
    clusters <- cutree(as.hclust(row_dend), k = 4)

##intersectBed -b data/cpgIslandExt.bed -a <( tail -n+2 results/2025-04-30/filt_anova.tsv | awk '$(NF)<0.05' ) |cut -f 1-4 > results/2025-04-30/incpgi.bed 
x=read.table("results/2025-04-30/incpgi.bed",header=F)
colnames(x)=c("chr","start","end","strand")
x$CpG=T
y=merge(df,x,all.x=T)
y$CpG[is.na(y$CpG)]=F
    right_anno=rowAnnotation(CpG= y[i,]$CpG)
    top_anno = HeatmapAnnotation(foo = anno_boxplot(m, height = unit(2, "cm")))
    png(file=paste0(o,"_heatmap.png"))
    Heatmap(m,top_annotation= top_anno,right_annotation=right_anno,row_split=clusters,show_row_names=F,cluster_columns=F,column_split=g)
    dev.off()

    for (k in unique(clusters)) {
      d <- df[i,][clusters == k, , drop = FALSE]
      write.table(d, file = paste0(o,"_cluster", k, ".tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
    }
'
}
pval=0.05
n=`tail -n+2 $odir/filt_anova.tsv | awk -v p=$pval '$(NF)<p' | wc -l`
N=`tail -n+2 $odir/filt_anova.tsv | wc -l `;
h=$odir/filt_heatmap.png;
echo '## Update : '`date`'
- Improved sensitivity in detecting '$n' / '$N' significant CpG sites (Pval <'$pval')
- Identified unbiased methylation patterns without applying prefiltering (previously required a 25% difference threshold)

## Results
- ANOVA results : [anova](/'$odir/filt_anova.tsv')
- Pval < '$pval' significant sites

| Heatmap | 
| :-: |
| ![hm](/'$h') |

| cluster_id | table |
| :-: | :-: |' > $o
for i in 1 2 3 5;do
    echo "| cluster_$i | [$i](/$odir/filt_cluster$i.tsv) |" >> $o
done

git add -A
git commit -am wooutlier
git push

exit
anno(){
i=filtered_anova.tsv
o=${i%.tsv}_anno
annotatePeaks.pl $i mm10 -go ${o}_go -annStats $o.stats > $o.tsv 
for i in 1 2 3 4;do
    continue;
    #i=merged_anova_pval005_cluster$i.tsv
    i=filtered_anova_cluster$i.tsv
    o=${i%.tsv}_anno
    annotatePeaks.pl $i mm10 -go ${o}_go -annStats $o.stats > $o.tsv 
done
}

bam2tsv(){
. ../../src/meth.sh
for f in *.bam;do
    n=${f%.dedu*}
    #bisbam2meth $f $n $n.tsv
    #chr1	4337982	4337982	+	11	8	3
done
}

cpg(){
    input=(
    2-year_1.tsv	2-year_3.tsv	E18pt5_1.tsv	E18pt5_3.tsv	Week4_1.tsv	Week4_3.tsv	
    2-year_2.tsv	2-year_4.tsv	E18pt5_2.tsv	E18pt5_4.tsv	Week4_2.tsv	Week4_4.tsv	
    )

    for f in ${input[@]};do
        n=${f%.tsv}
        #chr	start	end	strand	coverage	numCs	numTs
        tail -n+2 $f |  awk -v n=$n -v OFS="\t" '{print $1,$2,$2+1,n,$5,$6;}' |\
        intersectBed -a ../../data/cpgIslandExt.bed -b stdin -wa -wb |\
        sort -k1,1 -k2,3n | groupBy -g 1,2,3,8 -c 9,10, -o sum,sum
    done  |  awk '{print $1":"$2"-"$3"\t"$4"\t"($6/$5);}' | ca rc2mat -  > cpg.tsv
}

input=(
    E18pt5_1.tsv E18pt5_2.tsv E18pt5_3.tsv E18pt5_4.tsv
    Week4_1.tsv Week4_2.tsv Week4_3.tsv Week4_4.tsv
    2-year_1.tsv 2-year_2.tsv 2-year_3.tsv 2-year_4.tsv
)
. ../../src/util.sh;
#merge ${input[@]}  > merged.tsv

