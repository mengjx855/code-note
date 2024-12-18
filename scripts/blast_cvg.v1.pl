#!/usr/bin/perl
use warnings;
use strict;


die "perl $0 [in.blast] [out.blast] [%iden] [min.assign.len]\n" unless @ARGV == 4;
my ($in_f, $out_f, $iden, $minlen) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;
### print STDERR "Program $0 Start...\n";

my ($que, $sub, $sq, $ss, $nq, $ns) = ("","", 0, 0, 0, 0);
my (%que, %sub) = ();

open IN, $in_f or die $!;
open OT, ">$out_f" or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;
	
	next unless $s[2] >= $iden and $s[3] >= $minlen;

	($s[6], $s[7]) = ($s[7], $s[6]) unless $s[6] <= $s[7];
	($s[8], $s[9]) = ($s[9], $s[8]) unless $s[8] <= $s[9];
	
	my $similar = 1 - ($s[4] + $s[5]) / $s[3];

	if($que eq $s[0] and $sub eq $s[1]){
		for($s[6]..$s[7]){ next if exists $que{$_}; $que{$_}++; $nq++; $sq += $similar; }
		for($s[8]..$s[9]){ next if exists $sub{$_}; $sub{$_}++; $ns++; $ss += $similar; }
	}else{
		if($que ne ""){ printf OT "$que\t$sub\t$nq\t%.2f\t$ns\t%.2f\n", $sq/$nq*100, $ss/$ns*100;}

		($que, $sub, $sq, $ss, $nq, $ns) = ($s[0], $s[1], 0, 0, 0, 0);
		undef %que;
		undef %sub;
		for($s[6]..$s[7]){ $que{$_}++; $nq++; $sq += $similar; }
		for($s[8]..$s[9]){ $sub{$_}++; $ns++; $ss += $similar; }
	}
}
if($que ne ""){ printf OT "$que\t$sub\t$nq\t%.2f\t$ns\t%.2f\n", $sq/$nq*100, $ss/$ns*100;}
close IN;
close OT;

print STDERR "Program End...\n";
############################################################
sub function {
	return;
}
############################################################

