ref2bed12(){
    cat $1 | perl -ne 'chomp;
    my ($i,$v,$c,$t,$s,$e,$ts,$te,$n,$bs,$be,$dummy,$g)=split/\t/,$_; 
    $ts = $ts -1 ;
    my @ss=split/,/,$bs;
    my @ee=split/,/,$be;
    my @bsize = map { $ee[$_] - $ss[$_] } 0..($n-1);
    my @bstart = map { $ss[$_] - $s } 0..($n-1);
    print join("\t",$c,$s,$e,$v,0,$t, ### duplicated $g caused error 
    $ts,$te,"0,0,0",$n,join(",",@bsize),join(",",@bstart)),"\n";
'
}
