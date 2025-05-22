base=/mnt/vstor/SOM_GENE_BEG33/data/emseq
#ssh hxk728@pioneer.case.edu "find $base/bigdata/emseq-2025-04-11 -type d -name out" | sed "s|^$base/||" > out_dirs.txt
rsync -arv --relative --files-from=out_dirs.txt hxk728@pioneer.case.edu:$base/ bigdata/methylseq/
#scp hxk728@pioneer.case.edu:/mnt/vstor/SOM_GENE_BEG33/lscratch/hxk728/2-year_1/nf_work/2c/cdd559eaf36e047f70b136b2393c6b/BismarkIndex/*  bigdata/mm10/





