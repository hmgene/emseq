
cutn(){
local tmp=`mktemp -d`;
cat $1 > $tmp/a
Rscript <( echo '
    library(data.table)
    cols=strsplit("'$2'",",")[[1]];
    tt=fread("'$tmp/a'"); 
    fwrite(tt[,..cols],"'$tmp/b'",sep="\t",na = "NA",quote=F)
')
tail -n+2 $tmp/b 
}

merge(){
python <( echo 'import sys
import pandas as pd
fs=sys.argv[1:]
tt=None;
for f in fs:
    n=f.replace(".tsv","")
    d=pd.read_csv(f,sep="\t")
    d.columns = list(d.columns[:4]) + [f"{n}.{i}" for i in d.columns[4:]]
    if tt is None:
        tt = d
    else:
        tt = tt.merge(d, how="outer") 

tt = tt.fillna(0)
for col in tt.columns[4:]:
    tt[col] = tt[col].astype(int)

tt.to_csv(sys.stdout,sep="\t",index=False)
') $@
}
