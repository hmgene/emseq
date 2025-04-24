
for i in {1..5};do
#TermID	Term	Enrichment	logP	Genes in Term	Target Genes in Term	Fraction of Targets in Term	Total Target Genes	Total Genes	Entrez Gene IDs	Gene Symbols
#GO:0061646	positive regulation of glutamate neurotransmitter secretion in response to membrane depolarization	0.000343737107427677	-7.97563341535179	2	1	0.2	D
    cut -f 2,4,11 results/2025-04-11/*cluster$i*_go/biological_process.txt |\
    perl -ne 'chomp; my ($t,$p,$g)=split /\t/,$_;
    if(exp($p) < 0.005){ 
            $t=~s/regulation of/reg.of/g;
            $t=~s/positive reg.of/p.reg.of/g;
            $t=~s/negative reg.of/n.reg.of/g;

            my @tt=split/\s+/,$t;
            my @gg=split/,/,$g;
            map { 
                #print join("\t", join("_","cluster'$i'",@tt[0..2]), $_,1),"\n";
                print join("\t", join("_",,@tt[0..3]),"cluster'$i'", 1),"\n";
            } @gg;
    }
    ' 

done| ca rc2mat -  | Rscript -e 'tt=read.table("stdin",header=T,sep="\t")
library(ComplexHeatmap);
m=as.matrix(tt[,-1]);row.names(m)=tt[,1];
png("a_heatmap.png",height=1200)
Heatmap(m)
dev.off();
'
 
