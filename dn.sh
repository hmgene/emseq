

o=README.md
echo "### Table of Results" > $o
echo "|sample|summary|bedGraph|" >> $o
echo "|:--:|:--|:--|" >> $o
ssh hxk728@pioneer.case.edu ls -d '/mnt/vstor/SOM_GENE_BEG33/data/emseq/emseq/20250219_*/bismark' |\
while read d;do
    s=${d##*emseq/};s=${s%/bismark*};
    #scp hxk728@pioneer.case.edu:"$d/reports/*.html" data/${s}_reports.html
    #scp hxk728@pioneer.case.edu:"$d/methylation_calls/bedGraph/*.bedGraph.gz" data/${s}.bedGraph.gz 
    echo "| $s|[summary](https://raw.githack.com/hmgene/emseq/main/data/${s}_reports.html ) | [bedGraph](https://raw.githack.com/hmgene/emseq/main/data/${s}.bedGraph.gz ) |" >> $o
done 



