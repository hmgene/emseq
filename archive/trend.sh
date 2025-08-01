Rscript <( echo '
    t=read.table(text="
    Lin28b
    Hmga2
    Hic2
    Igf2bp2
    Igf2bp3
    ")
    tt=fread("results/filtered.3x.10bp.anova.anno.tsv")
    tt1=tt[`Gene Name` %in% t$V1] 
    p= as.matrix(tt1[,grep("_perc",names(tt1),value=T),with=F] )
    rownames(p) <- make.unique(paste0(tt1$`Gene Name`,".",substr(tt1$Annotation,1,10),".",tt1$chrom,":",tt1$start,"-",tt1$end))
    pdf("results/target_heatmap.pdf")
    Heatmap(p, row_names_gp = grid::gpar(fontsize = 5))
    dev.off()
    pdf("results/target_heatmap_scaled.pdf")
    Heatmap(t(scale(t(p))), row_names_gp = grid::gpar(fontsize = 5))
    dev.off()

    ###vln trend
    library(ggplot2)
    library(data.table)
    library(reshape2)

    tt[, `:=` (
      trend_EW = fifelse(mean_W > mean_E, "up", fifelse(mean_W < mean_E, "dn", "no_change")),
      trend_WY = fifelse(mean_Y > mean_W, "up", fifelse(mean_Y < mean_W, "dn", "no_change"))
    )]
    # Now classify the full trend E → W → Y
    tt[, trend_class := paste(trend_EW, trend_WY, sep = "_")]
    fwrite(tt,file="results/filtered.3x.10bp.anova.anno.trend.tsv",sep="\t")

    # Melt data to long format for ggplot (E, W, Y as condition)
    tt_long <- melt(tt, id.vars = c("group","trend_class","Gene Name"), measure.vars = c("mean_E", "mean_W", "mean_Y"),
                    variable.name = "condition", value.name = "methylation")
    # Clean up condition names
    tt_long[, condition := gsub("mean_", "", condition)]
    trend_counts <- tt[, .N, by = trend_class]
    trend_counts[, trend_label := paste0(trend_class, "\n(n=", N, ")")]

    tt_long <- merge(tt_long, trend_counts[, .(trend_class, trend_label)], by = "trend_class")
    group_counts <- tt_long[, .N, by = group]
    group_counts[, group_label := paste0(group, " (n=", N, ")")]
    group_label_map <- setNames(group_counts$group_label, group_counts$group)
    tt_long <- as.data.table(tt_long)
    tt_long[, group := group_label_map[group]]

    ggplot(tt_long, aes(x = trend_class, y = methylation, fill = condition)) +
      geom_violin(trim = FALSE, scale = "width", position = position_dodge(width = 0.8)) +
      stat_summary(fun = median, geom = "point", shape = 23, size = 2,
                   position = position_dodge(width = 0.8)) +
      facet_wrap(~ group, scales = "free_y", ncol = 1) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Methylation Trends by Group",
           x = "Trend (E→W→Y)", y = "Methylation (%)")


    ## summary
    dt=tt
    group_counts <- dt[, .N, by = group]
    dt <- merge(dt, group_counts, by = "group", suffixes = c("", "_group"))
    dt[, group := paste0(group, " (n=", N_group, ")")]
    dt[, N_group := NULL]

    trend_summary <- dt[, .(
      N = .N,
      mean_E = mean(mean_E, na.rm = TRUE),
      mean_W = mean(mean_W, na.rm = TRUE),
      mean_Y = mean(mean_Y, na.rm = TRUE)
    ), by = .(group )]

    setorder(trend_summary, group )
    print(trend_summary)
    #          group     N     mean_E     mean_W     mean_Y
    #         <char> <int>      <num>      <num>      <num>
    #1:         <NA>    14  0.2340333  0.2077105  0.5092456
    #2:           3'  1405 63.6838120 62.1841517 57.6932052
    #3:           5'   238 54.6511712 55.0127968 48.1433496
    #4:   Intergenic 57719 66.9876675 64.7996101 65.0313108
    #5:          TTS  1813 63.2491014 61.6091227 59.4552735
    #6:         exon  3712 55.9893814 62.6223668 65.0061654
    #7:       intron 50188 66.3645757 62.1072053 56.7808921
    #8:   non-coding   944 48.3941965 52.6688076 54.0259305
    #9: promoter-TSS  3910 53.7976471 54.1951853 53.0326458


    trend_counts <- dt[, .N, by = .(group, trend_class)]
    trend_wide <- dcast(trend_counts, group ~ trend_class, value.var = "N", fill = 0)
    print(trend_wide)
    #         group dn_dn dn_nc dn_up nc_dn nc_up up_dn up_nc up_up
    #1           3'   373     0   356     2     1   379     4   290
    #2           5'    66     1    52     0     0    68     1    50
    #3         exon   659     1   896     2     1  1027    12  1114
    #4   Intergenic 11706    40 18309    99    41 16211   211 11102
    #5       intron 14279    49 13243    95    25 13484   145  8868
    #6   non-coding   156     1   159     2     2   183     2   439
    #7 promoter-TSS   863     3  1072     2     2  1140    13   815
    #8          TTS   458     1   464     2     2   473     6   407
    #9         <NA>     1     0     7     0     0     1     0     5
    ')
}



#mkdir -p results/trend/
#rm results/trend/*
#cat  results/filtered.3x.10bp.anova.anno.trend.tsv | hm cutn - chrom,start,end,trend_class| cut -f 1-3,6 |\
#tail -n+2 |  awk -v OFS="\t" -v o="results/trend" '{out=o"/"$4; print $0 >> out}'

input=(
results/trend/dn_dn
results/trend/dn_nc
results/trend/dn_up
results/trend/nc_dn
results/trend/nc_up
results/trend/up_dn
results/trend/up_nc
results/trend/up_up
)
#parallel "annotatePeaks.pl {} mm10 -go {}_go -annStats {}.stats" ::: ${input[@]}
fn(){
    hm homer-go-sum  ${1}_go/biological_process.txt 0.0000001 40 50 |\
    perl -ne 'chomp; $_=~s/regulation of/reg./g; $_=~s/positive /p./g; $_=~s/negative /n./g; $_=~s/ /_/g;print join("\t","'${1##*/}'",$_,),"\n";'
};export -f fn
parallel --line-buffer fn {} ::: ${input[@]} | awk -v OFS="\t"  '{print $2,$1,1;}' | hm rc2mat - | perl -npe '$_=~s/\t(\d+)/\t1/g;$_=~s/nan/0/g;' | hm heatmap - o.pdf  7 14
exit









