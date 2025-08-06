odir="results/2025-08-06/";

summ(){
Rscript -e '
library(data.table)
idir="bigdata/2025-08-01/";
odir="results/2025-08-06/";
tt=fread(paste0(idir,"filtered_3rep_3x_per_group.csv.gz"))

CpG = colSums(tt[, grep("\\.CpG",names(tt),value=T), with=F],  na.rm = TRUE)
uCpG = colSums(tt[, grep(".uCpG",names(tt),value=T), with=F],  na.rm = TRUE)
id=substr(names(CpG),1,2);
group = substr(id,1,1);

df <- data.frame(id=id, CpG=CpG, uCpG=uCpG, group=group)
fwrite(df,paste0(odir,"/cpg_prop_by_group.csv"))
'
}

cd $odir
Rscript -e 'library(rmarkdown);render("README.Rmd")'
cd -
git commit -am 1 
git push
