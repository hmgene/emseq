bismark-anova () 
{ 
    usage="$FUNCNAME <.cov[.gz]> <meta.tsv>";
    if [ $# -lt 2 ]; then
        echo "$usage";
        return;
    fi;
    local tmp=`mktemp -d`;
    Rscript -e 'i="'$1'"; m="'$2'";
    library(data.table);
    #i="bigdata/filtered.cov.gz"; m="meta.txt";
    tt=fread(i); me=fread(m)
    p = tt[, grep("_perc",names(tt),value=T),with=F]
    names(p)=gsub("_perc","",names(p))
    g=me[match(names(p),me[,sample]),group]
    res = cbind(tt, t( apply(p, 1, function(x){
      r <- try(aov(x ~ g, na.action = na.omit), silent = TRUE)
      if (inherits(r, "try-error")) return(c(pval = NA, gW = NA, gY = NA))
      pval <- summary(r)[[1]]$`Pr(>F)`[1]
      coefs <- coef(r)
      gW <- if ("gW" %in% names(coefs)) coefs["gW"] else 0
      gY <- if ("gY" %in% names(coefs)) coefs["gY"] else 0
      return(c(pval = pval, gW, gY))
    })));
    fwrite(res,file="'$tmp/a'",sep="\t")
'
}
