run-methylseq(){
export PATH=/mnt/vstor/SOM_GENE_BEG33/mamba/miniforge3/bin:$PATH
input="bigdata/emseq-2025-04-11/*/*_R1.fastq.gz" ## need to be merged

for f in $input;do
        s=${f##*/};s=${s%_R1.fastq*}
        odir=${f%/*}
cd $odir
echo "#!/bin/bash
module load singularity
nextflow run nf-core/methylseq -r 2.5.0 --genome mm10 --input input.csv --outdir ./out --em_seq -profile singularity \
-w /mnt/vstor/SOM_GENE_BEG33/lscratch/hxk728/$s/nf_work/
" | sbatch --mem=128g -c 24 --time 100:00:00
cd -
}



