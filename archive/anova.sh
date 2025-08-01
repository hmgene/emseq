input=(
filtered.3x.cov.gz  filtered.5x.cov.gz  
)

for i in ${input[@]};do
for b in 1 10;do
        echo "#!/bin/bash
        hm bismark-cov-bin $i $b | hm bismark-anova - meta.tsv 0.05 | gzip -c > ${i%.cov.gz}.${b}bp.anova.tsv.gz
        " #| sbatch --mem=64g -c 8
done
done

