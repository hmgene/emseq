Rscript -e '
library(data.table)
idir="bigdata/2025-08-01/";
odir="results/2025-08-06/";
tt=fread(paste0(odir,"filtered_3rep_3x_per_group.csv.gz"))

CpG = colSums(tt[, grep(".CpG",names(tt),value=T), with=F],  na.rm = TRUE)
uCpG = colSums(tt[, grep(".uCpG",names(tt),value=T), with=F],  na.rm = TRUE)
group = factor(substr(names(CpG),1,1))

df <- data.frame(CpG=CpG, uCpG=uCpG, group=group)
fwrite(df,paste0(odir,"/cpg_prop_by_group.csv"))



'

