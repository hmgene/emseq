input=(
bigdata/meth/10_2-year_2.tsv  bigdata/meth/12_2-year_4.tsv  bigdata/meth/2_E18pt5_2.tsv  bigdata/meth/4_E18pt5_4.tsv  bigdata/meth/6_Week4_2.tsv  bigdata/meth/8_Week4_4.tsv
bigdata/meth/11_2-year_3.tsv  bigdata/meth/1_E18pt5_1.tsv   bigdata/meth/3_E18pt5_3.tsv  bigdata/meth/5_Week4_1.tsv   bigdata/meth/7_Week4_3.tsv  bigdata/meth/9_2-year_1.tsv
)

#chr	start	end	strand	coverage	numCs	numTs
#chr1	3057454	3057454	-	10	10	0
#chr1	3135943	3135943	-	10	10	0
for f in ${input[@]};do
    n=${f##*/};n=${n%.tsv};
    tail -n+2 $f | awk -v n=$n '{print $1":"$2$4"\t"n"\t"($6/$5);}'  
done | ca rc2mat - > bigdata/meth/merged_perc.tsv


R --no-save -e '
tt=read.table("bigdata/meth/merged_perc.tsv",header=T);
m=tt[,-1];row.names(m)=tt[,1]
g=rep("E",ncol(m)); g[grep("Week",colnames(m))] = "W"; g[grep("year",colnames(m))] = "Y"

v=apply(m[,g=="W"],1,median)
m1=m[order(-v)[1:1000],]
library(ComplexHeatmap)

Heatmap(m1,column_split=g)


i=apply(m,1,mean)

colnames(m) <- paste0(colnames(m), "_", g)
m_df <- as.data.frame(m[i>0.25,])
m_df$ID <- rownames(m_df)
m_long <- melt(m_df, id.vars = "ID", variable.name = "Sample", value.name = "Methylation")
m_long$Group <- sub(".*_(\\w+)$", "\\1", m_long$Sample)

# Plot
ggplot(m_long, aes(x = Sample, y = Methylation, fill = Group)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "Violin Plot per Sample (split, colored by group)",
       x = "Sample", y = "Methylation Level") +
  scale_fill_brewer(palette = "Set2")
'


exit


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



