#!/usr/bin/perl
use warnings;
use strict;

############################################################
### Author: lishenghui@genomics.org.cn                   ###
### Version: 1.0 (Apr. 1st, 2012)                        ###
############################################################

die "perl $0 [genelist] [reference.fa] [genelist.fa (outfile)]\n" unless @ARGV == 3;
my ($in_f, $ref_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;
### print STDERR "Program $0 Start...\n";

my %gene;


if ($in_f =~ /\.gz$/) { open IN, "gzip -dc $in_f |" or die $!;
} else { open IN, $in_f or die $!;
}
while(<IN>){
	chomp;
	s/\s+.*$//g;
	$gene{$_} = $1;
}
close IN;

my $bool = 0;

if ($ref_f =~ /\.gz$/) { open IN, "gzip -dc $ref_f |" or die $!;
} else { open IN, $ref_f or die $!;
}
open OT, ">$out_f" or die $!;
while(<IN>){
	chomp;
	next if $_ eq "";
	if(/^>(\S+)/){
		if(exists $gene{$1}){ $bool = 1; delete $gene{$1};
		} else { $bool = 0;
		}
	}
	print OT "$_\n" if $bool == 1;
}
close IN;
close OT;

print STDERR "Program End...\n";
############################################################
sub function {
	return;
}
############################################################

