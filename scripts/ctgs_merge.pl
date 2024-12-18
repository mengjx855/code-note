#!/usr/bin/perl
use warnings;
use strict;

die "perl $0 <ctgs_list> <out_f>" unless @ARGV eq 2;
my ($in_f, $out_f) = ($ARGV[0], $ARGV[1]);

open IN, "<$in_f" or die "Can't open $in_f: $!\n";
open OU, ">$out_f" or die "Can't open $out_f: $!\n";

while (<IN>) {
    
    chomp;
    open I, "<$_";
    my @s = split/\//;
    $s[-1] =~ s/.contigs.fa//;
    
    my $i = 1;
    while (<I>) {
        chomp;
        
        if (/>/) {
            print OU ">$s[-3]_$s[-1]_$i\n";
            $i++;
        } else { print OU "$_\n"; }   
    }
}
close IN;
close OU;