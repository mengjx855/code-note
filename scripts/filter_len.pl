#!/usr/bin/perl
use warnings;
use strict;

die "usage: perl $0 [old.fa] [new.fa] [length]\n" unless @ARGV == 3;

my ($in_f, $out_f, $length) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

my $break = 100;

if($in_f =~ /\.gz$/) { open IN, "gzip -dc $in_f |" or die $!;
} else { open IN, $in_f or die $!;
}
open OT, ">$out_f" or die $!;

my ($seq, $head) = ("", "");
while (<IN>) {
    chomp;

    if (/^(>\S+)/) {
		&print_seq($seq) unless $seq eq "";
		($seq, $head) = ("", $1);
		next;
	}
	s/n/N/g; #########
	$seq .= $_;
	
}
&print_seq($seq) unless $seq eq "";
close IN;
close OT;
############################################################
sub print_seq {
	my $seq = shift;
	my $seq_len = length $seq;

	my $seu = $seq; 
	$seu =~ s/N+//g;
	my $seu_len = length $seu;

	return unless $seu_len >= $length; #########
	
	my $line = "";
	for (my $i = 0; $i < $seq_len; $i += $break) {
		$line .= substr($seq, $i, $break)."\n";
	}
	chomp $line;

	### printf OT "$head length=%d/%d\n$line\n", $seu_len, $seq_len;
	print OT "$head\n$line\n";
}
############################################################

