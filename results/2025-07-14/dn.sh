i="
BMHSC_RNA-Seq_Rep1_GSE122908  SRR8241576
BMHSC_RNA-Seq_Rep2_GSE122908  SRR8241577
FL_HSC_rep1_GSE104689 SRR6144580
FL_HSC_rep2_GSE104689 SRR6144581
"
odir=/mnt/vstor/SOM_GENE_BEG33/data/emseq/fastq
echo "$i" | grep -v "^$" | grep -v "^#" |  while read -r x y;do
mkdir -p $odir
echo "#!/bin/bash
	cd $odir
	hm fastq-dump --gzip  --split-3 $y
	cd -
" | sbatch
done 
