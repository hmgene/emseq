
anova(){
usage="$FUNCNAME <input.tsv> <bin=1> <strand=1>" 
if [ $# -lt 3 ];then echo "$usage";return;fi
input=$1; bin=$2; strand=$3; output=$4

Rscript <( echo '
    input="'$input'"; o="'$output'"; bin_size = '$bin'; strand = '$strand';
    #input="bigdata/2025-04-11/merged.tsv"; o="results/2025-05-08/cmethy_50bp.tsv"; bin_size=1; strand=1
    by=c("chr","bin"); if(strand){ by=c(by,"strand");}
    library(data.table)
    tt=fread(input); tt[, bin := floor(start / bin_size) ]; tt[is.na(tt)]=0;
    nC <- grep("numCs", names(tt), value = TRUE); nN <- grep("coverage", names(tt), value = TRUE)
    d = as.data.table(tt[, lapply(.SD, sum), by = by, .SDcols = c(nC,nN)] )
    d[, `:=` ( start=bin*bin_size, end=(bin+1)*bin_size )]

    for ( j in sub(".coverage","",names(d)[grep(".coverage",names(d))]) ){
        d[, (paste0(j,"_perc")) := 100 * get( paste0(j,".numCs") )/ get ( paste0(j,".coverage")) ]
    }
    
    ## filtering
    x=d[,rowSums(is.na(.SD)),.SDcols = grep("^E.+_perc",names(d),value=T)]; d=d[x<3,]
    x=d[,rowSums(is.na(.SD)),.SDcols = grep("^W.+_perc",names(d),value=T)]; d=d[x<3,]
    x=d[,rowSums(is.na(.SD)),.SDcols = grep("^2.+_perc",names(d),value=T)]; d=d[x<3,]
    
    p = d[, grep("_perc",names(d),value=T),with=F]
    ## ANOVA
    g=rep("E",ncol(p)); g[grep("^W",colnames(p))] = "W"; g[grep("^2",colnames(p))] = "Y"
    res = cbind(d, t( apply(p, 1, function(x){
      r <- try(aov(x ~ g, na.action = na.omit), silent = TRUE)
      if (inherits(r, "try-error")) return(c(pval = NA, gW = NA, gY = NA))
      pval <- summary(r)[[1]]$`Pr(>F)`[1]
      coefs <- coef(r)
      gW <- if ("gW" %in% names(coefs)) coefs["gW"] else 0
      gY <- if ("gY" %in% names(coefs)) coefs["gY"] else 0
      return(c(pval = pval, gW, gY))
    })))
    write.table(res,stdout(),sep="\t",col.names=T,row.names=F,quote=F)
')
}

multibin(){
usage="$FUNCNAME <output> <1bp> <10bp> <50bp>"
if [ $# -lt 3 ];then echo "$usage";return;fi;
local input=${@:2};
local output=$1;
#input=( results/2025-05-08/anova_1bp.tsv results/2025-05-08/anova_10bp.tsv results/2025-05-08/anova_50bp.tsv )

echo ${input[@]};return;
    Rscript <( echo 'input <- commandArgs(trailingOnly = TRUE);
    ## order by bin _\\d+bp
    output="'$output'" # output="results/2025-05-08/anova_multi"
    library(data.table)

    cn = c("chr","strand","start","end") ## last two 
    cn1 = c("pval","gW","gY","perc")
    tt=NULL;
    for( i in input){
        bz=sub(".*_(\\d+)bp.*", "\\1", i) 
        tmp= fread(i); j=c(cn,grep(paste(cn1,collapse="|"),names(tmp),value=T)); tmp=tmp[,..j]; tmp[, start := as.integer(start)]; tmp[, end := as.integer(end)]
        setkeyv(tmp, cn)
        setnames(tmp,c(cn,paste0(j[(length(cn)+1):(length(j))],"_",bz)))
        if( is.null(tt) ){ tt=tmp;
        }else{ 
            tt= foverlaps(tt,tmp,type="any",nomatch=NA); 
            tt[, start := NULL]; tt[, end := NULL]
            setnames(tt, gsub("^i\\.", "", names(tt)))
        }
    }
    fwrite(tt,paste0(output,".tsv"),sep="\t")

    pval_cols <- grep("^pval_", names(tt), value = TRUE)
    tt[, min_pval := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = pval_cols]
    m=as.matrix(tt[min_pval < 0.05, grep("perc",names(tt),value=T),with=F])
    m=t(scale(t(m)));
    m[is.na(m)]=0;
    gbin=sub(".+_(\\d+)$","\\1",colnames(m))
    gsam=sub("^([^_]+).+","\\1",colnames(m))

    library(ComplexHeatmap)
    png(file=paste0(output,"_heatmap.png"))
    Heatmap(m,use_raster=F,row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),column_split=gsam)
    dev.off();

    png(paste0(output,"_corbysample.png"))
    Heatmap(cor(m),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),column_split=gsam,row_split=gsam)
    dev.off();
    png(paste0(output,"_corbybin.png"))
    Heatmap(cor(m),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),column_split=gbin,row_split=gbin)
    dev.off();

    ## strandness
    y=tt[min_pval<0.05 , if (uniqueN(strand) > 1) .SD, by = .(chr, start, end)]
    y1=dcast(y, chr + start + end ~ strand, value.var = c("gW_1","gW_10","gW_50","gY_1","gY_10","gY_50"), fun.aggregate = mean)
    Heatmap(y1[,4:ncol(y1)],column_split=sub(".+_(\\d+)_[+-]","\\1",colnames(y1)[4:ncol(y1)]))
    ') ${@:2}
}
