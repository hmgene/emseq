odir="bigdata"; mkdir -p $odir

make-bg(){
    Rscript <( echo '
        library(data.table);
        input="../../bigdata/2025-08-01/filtered_3rep_3x_per_group.csv.gz"
        odir="bigdata"
        input_sig="sig.csv";
        tt=fread(input)
        tt[,`:=`(start = as.integer(start), end = as.integer(start+1))]
        tt[, grep("Y1|E2", names(tt),value=T) := NULL ]
        for (g in c("E", "W", "Y")) {
          tt[, paste0("t", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\.uCpG"), names(tt), value = TRUE)]
          tt[, paste0("c", g) := rowSums(.SD, na.rm = TRUE), .SDcols = grep(paste0(g, "\\d+\\.CpG"), names(tt), value = TRUE)]
          tt[, paste0("p", g) := get(paste0("c", g)) / (get(paste0("t", g))+get(paste0("c",g)))]
        }

        fwrite(tt[,.(chrom,start,end,cE)],paste0(odir,"/EM_cE.bedGraph.gz"),sep="\t",col.names=F)
        fwrite(tt[,.(chrom,start,end,tE)],paste0(odir,"/EM_tE.bedGraph.gz"),sep="\t",col.names=F)
        fwrite(tt[,.(chrom,start,end,cW)],paste0(odir,"/EM_cW.bedGraph.gz"),sep="\t",col.names=F)
        fwrite(tt[,.(chrom,start,end,tW)],paste0(odir,"/EM_tW.bedGraph.gz"),sep="\t",col.names=F)
        fwrite(tt[,.(chrom,start,end,pE)],paste0(odir,"/EM_pE.bedGraph.gz"),sep="\t",col.names=F)
        fwrite(tt[,.(chrom,start,end,pW)],paste0(odir,"/EM_pW.bedGraph.gz"),sep="\t",col.names=F)
    ')
}
#make-bg
fn(){
echo "
    gunzip -dc $1 | hm bg2bw - mm10 > ${1%.gz}.bw
"
};export -f fn

parallel fn {} ::: bigdata/EM_p{W,E}.bedGraph.gz
exit
#bigdata/EM_cE.bedGraph.gz  bigdata/EM_cW.bedGraph.gz  bigdata/EM_tE.bedGraph.gz  bigdata/EM_tW.bedGraph.gz 
#hm ucsc-refflat mm10 | hm ucsc-refflat2bed12 - 0  | hm bed5p -  | sort -u  > $odir/tss.bed

#awk '$4 > 1 && $8 < 0.05{ print $1"\t"$1;}'  Beaudin.txt   |  hm repl bigdata/tss.bed - 1 > bigdata/tss_FL_Beaudin.bed
#awk '$4 < -1 && $8 < 0.05{ print $1"\t"$1;}'  Beaudin.txt   |  hm repl bigdata/tss.bed - 1 > bigdata/tss_BM_Beaudin.bed

#bw=( bigdata/EM_cE.bedGraph.bw	bigdata/EM_cW.bedGraph.bw	bigdata/EM_tE.bedGraph.bw	bigdata/EM_tW.bedGraph.bw )
bw=( bigdata/EM_pE.bedGraph.bw	bigdata/EM_cW.bedGraph.bw	bigdata/EM_tE.bedGraph.bw	bigdata/EM_tW.bedGraph.bw )
bd=( bigdata/tss_FL_Beaudin.bed bigdata/tss_BM_Beaudin.bed  ) 
#computeMatrix reference-point -S ${bw[@]}  -R ${bd[@]} -a 3000 -b 3000 -o bigdata/EM_matrix -bs 50 
#plotHeatmap -m bigdata/EM_matrix -o EM_heatmap.png
