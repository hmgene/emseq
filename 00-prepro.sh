#for f in  bigdata/methylseq/*/out/bismark/methylation_calls/methylation_coverage/*deduplicated.bismark.cov.gz;do
#    s=`echo $f| cut -d "/" -f 3`;
#    l=${s:0:1}; l=${l/2/Y}; r=${s: -1}; s=$l$r
#    echo "$s $f"
#done 
Rscript -e '
odir="bigdata/2025-08-01/";
input="Y1 bigdata/methylseq/2-year_1/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y2 bigdata/methylseq/2-year_2/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y3 bigdata/methylseq/2-year_3/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
Y4 bigdata/methylseq/2-year_4/out/bismark/methylation_calls/methylation_coverage/2-year_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E1 bigdata/methylseq/E18pt5_1/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E2 bigdata/methylseq/E18pt5_2/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E3 bigdata/methylseq/E18pt5_3/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
E4 bigdata/methylseq/E18pt5_4/out/bismark/methylation_calls/methylation_coverage/E18pt5_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W1 bigdata/methylseq/Week4_1/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W2 bigdata/methylseq/Week4_2/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W3 bigdata/methylseq/Week4_3/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
W4 bigdata/methylseq/Week4_4/out/bismark/methylation_calls/methylation_coverage/Week4_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
input=read.table(text=input);
library(data.table)
min_cov=3
tt=NULL
for( i in 1:nrow(input)){
    tmp=fread(input[i,2],col.names=c("chrom","start","end","perc","CpG","uCpG"));
    tmp=tmp[ CpG + uCpG >= min_cov,.(chrom,start,CpG,uCpG) ];
    setnames(tmp, c("chrom","start", paste0(input[i,1],c(".CpG",".uCpG"))))

    j=setdiff(names(tmp),c("tot","end"))
    if(is.null(tt)){
        tt=tmp[,..j];
    }else{
        tt=merge(tt,tmp[,..j],all=T)
    }
}
fwrite(file=paste0(odir,"merged.csv.gz"),tt)

groups <- c("E", "Y", "W")
for (g in groups) {
  pattern <- paste0("^", g, "\\d+\\.(?:CpG|uCpG)$")  # e.g., Y1.CpG, E3.uCpG
  cols <- grep(pattern, names(tt), value = TRUE, perl = TRUE)
  if (length(cols) == 0) {
    tt[, paste0(g, "_nonNA") := 0]
  } else {
    tt[, paste0(g, "_nonNA") := rowSums(!is.na(.SD)), .SDcols = cols]
  }
}

min_g=3 * 2
filtered <- tt[
  (E_nonNA >= min_g) & (Y_nonNA >= min_g) |
  (E_nonNA >= min_g) & (W_nonNA >= min_g) |
  (W_nonNA >= min_g) & (Y_nonNA >= min_g)
]
filtered[, c(paste0(groups, "_nonNA")) := NULL]
fwrite(filtered, file=paste0(odir,"filtered_3rep_3x_per_group.csv.gz"))



'

