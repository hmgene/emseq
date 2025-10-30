
#samtools view ../../bigdata/E18pt5_1/bismark/deduplicated/E18pt5_1.deduplicated.sorted.bam |\
#head -n 10000 |\
#hm bismark-sam21 - 

metilen-pp(){
Rscript -e '
library(data.table)
tt=fread("/mnt/vstor/SOM_GENE_BEG33/emseq/2025-09-11/filtered_3rep_3x_per_group.csv")
for ( s in unique(gsub("^([EWY]\\d)\\..*$", "\\1", grep("^[EWY]\\d\\.", names(tt), value = TRUE), perl = TRUE)) ){
  tt[, paste0(substr(s,0,1),"_",s, ".rate") := get(paste0(s, ".CpG")) / (get(paste0(s, ".CpG")) + get(paste0(s, ".uCpG")))]
}
d=tt[, .SD, .SDcols = c("chrom", "start", grep(".rate", names(tt), value = TRUE))]
d[is.na(d)] <- "-"
setnames(d, "start", "pos")
fwrite(d,"../../bigdata/metilen/metilen_input.tsv",sep="\t")
'
};export -f metilen-pp

metilen-rn(){
i=../../bigdata/metilen/metilen_input.tsv
echo "#!/bin/bash
#metilen-pp
../../tools/metilene_v0.2-9/metilene -t 8 -m 2 -a W -b E $i > metilen_output_WvsE.tsv
../../tools/metilene_v0.2-9/metilene -t 8 -m 2 -a Y -b W $i > metilen_output_YvsW.tsv
" | sbatch  --mem=64g -c 16
}

merge-dmr(){
{
awk '$8<0.05' metilen_output_YvsW.tsv 
awk '$8<0.05' metilen_output_WvsE.tsv 
} | sort -k1,1 -k2,3n | mergeBed -i stdin > merged.bed
}

Rscript -e '
library(data.table)
m=fread("merged.bed")
setnames(m,c("chrom","start","end"))
tt=fread("../../bigdata/metilen/metilen_input.tsv")
setnames(tt, "pos", "start")  # rename for foverlaps
tt[, end:= start]  
setkey(tt, chrom, start,end)
setkey(m, chrom, start, end)

rate_cols <- grep("\\.rate$", names(tt), value = TRUE)
tt[, (rate_cols) := lapply(.SD, as.numeric), .SDcols = rate_cols] 
res <- foverlaps(tt, m, by.x = c("chrom", "start", "end"),
                       by.y = c("chrom", "start", "end"),
                       type = "within", nomatch = 0L)

agg <- res[, lapply(.SD, mean, na.rm = TRUE), by = .(chrom, start, end), .SDcols = rate_cols]
x=agg[, ..rate_cols]

library(ComplexHeatmap)
pdf("heatmap.pdf")
Heatmap(t(scale(t(x))))
dev.off();

pdf("heatmap_raw.pdf")
Heatmap(x)
dev.off();

'
