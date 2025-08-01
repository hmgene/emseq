
tail -n+2 results/summary_1k_diffmeth.tsv | cut -f 1 | sed "s/[\.-]/\t/g" |\
sort -k1,1 -k2,3n |  mergeBed -i stdin -d 10000 > results/target_1k.bed
