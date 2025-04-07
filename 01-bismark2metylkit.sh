
fn(){
    mkdir -p ${3%/*} 
    . src/meth.sh; bisbam2meth $1 $2 $3
};export -f fn;
#parallel -j 4 fn {} {/.} bigdata/meth/from_bisbam/{/.}.tsv ::: bigdata/bams/*.bam

for f in bigdata/meth/from_bisbam/*.tsv;do
    #bigdata/meth/from_bisbam/20250219_10_2-year_2_MY12882_S10_L006.deduplicated.sorted.tsv
    n=${f##*20250219_}; n=${n%_MY12*};
    echo $n;
done | sort -u  | while read -r n;do
    o=bigdata/meth/$n.tsv
    python <( echo '
import sys;
import pandas as pd;
input_files = sys.argv[1:]
dfs = []
for f in input_files:
    df = pd.read_csv(f, sep="\t")
    dfs.append(df)
combined = pd.concat(dfs, ignore_index=True)
merged = combined.groupby(["chr", "start", "end", "strand"], as_index=False).sum(numeric_only=True)
merged = merged.sort_values(by=["chr", "start"])
merged.to_csv(sys.stdout, sep="\t", index=False)
') bigdata/meth/from_bisbam/*$n* > $o
done

