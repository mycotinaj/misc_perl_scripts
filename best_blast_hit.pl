#! /usr/bin/perl -w

#prints single best blast hit (highest bit score) for each of the query sequences in a blast output (tabular format)
# USAGE:  cat <blastfile> | best_blast_hit.pl > <outfile> 


while(<>){
	chomp;
	@spl=split(/\t/, $_);
	$bit{$spl[0]} = 0 unless defined($bit{$spl[0]});
	$best{$spl[0]}=$_ unless $bit{$spl[0]} > $spl[11]; 
	$bit{$spl[0]}=$spl[11] unless $bit{$spl[0]}>$spl[11];
}
foreach $key (keys %best){
	print "$best{$key}\n";
}
