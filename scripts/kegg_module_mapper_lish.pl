#!/usr/bin/perl
use warnings;
use strict;
die "Usage: perl $0 [ko.list] [ko.list.out]\n" unless @ARGV == 2;
my ($in_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

my $mod_f = "/share/data1/database/KEGG/KEGG20230401/module_map.reaction";

my %ko = ();
open IN, $in_f or die $!;
while(<IN>){
	chomp;
	if(/^(K\d\d\d\d\d)/){ $ko{$1}++; }
}
close IN;
my @k = keys %ko; print STDERR "Number-of-KOs: ".($#k+1)."\n";

my (%best, %best_v) = ();
print STDERR "read file $mod_f\n";
open IN, $mod_f or die $!;
while(<IN>){
	chomp;
	my ($mod,$path,$list) = split /\t+/;
	my @list = split / /,$list;
	my $sum = 0;
	
	for my $i(0..$#list){ 
		my $k = $list[$i];

		if($k=~/^K\d\d\d\d\d$/){	
			next unless exists $ko{$list[$i]};
			$sum++;
			$list[$i] .= "(m)";
		}else{
			$k =~ s/-K\d\d\d\d\d//g;
			my @new = split /\+/, $k;
			my $new = 1;
			for(@new){ $new=0 unless exists $ko{$_}; }
			next unless $new ==1;
			$sum++;
			$list[$i] = "".(join "+",@new)."(m)";
		}
	}	
	my $pct = int($sum/@list*10000+0.5)/100;

	unless(exists $best_v{$mod} and $best_v{$mod} >= $pct){
		###$best{$mod} = "$path\t$sum\/".($#list+1)."\t".(join " ",@list);
		$best{$mod} = "$sum\/".($#list+1)."\t".(join " ",@list);
		$best_v{$mod} = $pct;
	}
}
close IN;

open OT, ">$out_f" or die $!;
for(sort keys %best){
	print OT "$_\t$best_v{$_}\t$best{$_}\n";
}
close OT;
print STDERR "Program End...\n";
############################################################

