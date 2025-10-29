i="
BMHSC_RNA-Seq_Rep1_GSE122908  SRR8241576
BMHSC_RNA-Seq_Rep2_GSE122908  SRR8241577
FL_HSC_rep1_GSE104689 SRR6144580
FL_HSC_rep2_GSE104689 SRR6144581
"
r=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/mm10/scripts/rnaseq_mm10_v1.sh
#r=`realpath ./rnaseq_mm10_v1.sh`

odir=/mnt/vstor/SOM_GENE_BEG33/data/emseq/fastq
ddir=/mnt/vstor/SOM_GENE_BEG33/RNA_seq/mm10/DATA/

ls -la $odir;
echo "$i" | grep -v "^$" | grep -v "^#" |  while read -r x y;do
mkdir -p $odir
	l1=$odir/${y}_1.fastq.gz; l2=${y}_2.fastq.gz;
	r1=${x}_R1.fastq.gz; r2=${x}_R2.fastq.gz;
echo "#!/bin/bash
	cd $odir
	#hm fastq-dump --gzip  --split-3 $y
	if [ ! -f $r1 ];then mv $l1 $r1; fi
	if [ ! -f $r2 ];then mv $l2 $r2; fi
	if [ -f $r1 -a -f $r2 ];then
		$r $x $r1 $r2
	fi
	cd -
" #| sbatch --mem=24g -c 16
done 
