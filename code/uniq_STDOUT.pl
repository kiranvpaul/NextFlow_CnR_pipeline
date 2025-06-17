#!/usr/bin/perl
BEGIN { push @INC, 'lib'}

#Author: Srinivas Ramachandran

use genome_size;
use ngs;

die "Usage: perl uniq.pl <BED FILE> > OUT_BED\n" if(!$ARGV[0]);

if($ARGV[0]=~/gz$/){
	open(FILE,"gunzip -c $ARGV[0] |") || die "cannot open pipe to $ARGV[0] $! \n";
}else{
	open(FILE,$ARGV[0]) || die "cannot open $ARGV[0] $!\n";
}

while(chomp($line=<FILE>)){
	$lno++;
	@temp = split /[\ \s\n\t]+/, $line;
	$frag_key = $temp[0]."_".$temp[1]."_".$temp[2]."_".$temp[5];
	$reads{$frag_key} = $line;
}

@temp = keys(%reads);

$n_uniq = $#temp+1;

$pct_uniq = int($n_uniq*1000/$lno)/10;

print STDERR "$ARGV[0] Total_Reads: $lno Unique_reads: $n_uniq percent_uniq: $pct_uniq\n";

foreach $i (keys(%reads)){
	print "$reads{$i}\n";
}
