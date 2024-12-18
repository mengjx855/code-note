#!/share/software/perl/bin/perl
use warnings;
use strict;

############################################################
### Author: lishenghui@genomics.org.cn                   ###
### Version: 1.0 (Apr. 1st, 2012)                        ###
############################################################

die "perl $0 [blastout (input)] [outfile]\n" unless @ARGV == 2;
my ($in_f, $out_f) = @ARGV;
### print STDERR "Program $0 Start...\n";

my $acc2taxid_f = "/home/zhangy/bin/.depends/ncbi/nucl_gb.accession2taxid.gz";
my $node_f = "/home/zhangy/bin/.depends/ncbi/nodes.dmp";
my $name_f = "/home/zhangy/bin/.depends/ncbi/names.dmp";

my %taxa = ();
$taxa{"superkingdom"} = "k";
$taxa{"phylum"} = "p";
$taxa{"class"} = "c";
$taxa{"order"} = "o";
$taxa{"family"} = "f";
$taxa{"genus"} = "g";
$taxa{"species"} = "s";

############################################################

my (%node, %rank, %name, %tax) = ();

open IN, $node_f or die $!;
while(<IN>){
	chomp;
	my @s = split /\t+/;
	$node{$s[0]} = $s[2];
	$rank{$s[0]} = $s[4];
}
close IN;

open IN, $name_f or die $!;
while(<IN>){
	chomp;
	next unless /scientific name/;

	my @s = split /\t+/;
	$name{$s[0]} = $s[2];
}
close IN;

open IN, "gzip -dc $acc2taxid_f |" or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;
	$tax{$s[1]} = $s[2];
}
close IN;

print STDERR "processing...\n";
############################################################

my (%gene, %ident, %species) = ();

open IN, $in_f or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;

	my $t = $tax{$s[1]}; 
	next if $t eq "notax"; #########
	
	push @{$gene{$s[0]}}, $t;
	push @{$ident{$s[0]}}, $s[2];
	$species{$t} = 1;
}
close IN;

my (%s_node, %s_rank, %s_name, %info) = ();

foreach my $s(keys %species){

	my (@node, @rank, @name) = ();

	my $ss = $s;
	while(exists $node{$ss} and $ss != $node{$ss}){ ######
		my $r = $rank{$ss};
		if(exists $taxa{$r}){
			unshift @node, $ss;
			unshift @rank, $rank{$ss};
			unshift @name, $taxa{$r}."_".$name{$ss};
		}
		$ss = $node{$ss};
	}

	###
	if(@node > 0){
		$s_node{$s} = "".(join ";", @node);
		$s_rank{$s} = "".(join ";", @rank);
		$s_name{$s} = "".(join ";", @name);
		$info{$s_node{$s}} = $s_name{$s};
	}
}

print STDERR "writing...\n";
############################################################

open OT, ">$out_f" or die $!;
foreach my $g(sort keys %gene){

	my @tax = @{$gene{$g}};
	my @idt = @{$ident{$g}};

	my %spe = ();
	for(@tax){
		$spe{$s_node{$_}}++;
	}
	my @spe = sort{$spe{$b}<=>$spe{$a}} keys %spe;
		
	my $s = $spe[0];

	my $max = 0;
	for(0..$#tax){
		next unless $s_node{$tax[$_]} eq $s;
		$max = $idt[$_] if $idt[$_] > $max;
	}
	
	my $sum = $#tax + 1;
	my $pct = int ($spe{$s}/$sum * 1000 + 0.5)/10;

	print OT "$g\t$spe{$s}/$sum\t$pct\t$max\t$s\t$info{$s}\n";
}
close OT;

print STDERR "Program End...\n";
############################################################
sub min {
	return $_[0] if $_[0] <= $_[1];
	return $_[1];
}
############################################################

