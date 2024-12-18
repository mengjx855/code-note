#!/usr/bin/perl
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-09-27, 21:17:40
# modified date: 2023-11-21, 21:17:40
use warnings;
use strict;

my $db = "/share/data1/database/KEGG/KEGG20230401/KEGG_API/KO2Pathway.tsv";
my (%map, %count, %sum) = ();

open DB, "<$db";
while (<DB>) {
    chomp;
    my @s = split/\t/;
    push @{$map{$s[1]}}, $s[0];
}

for my $i (keys %map){
    $count{$i} = 0;
}

while (<>) {
    chomp;
    for my $k (keys %map){
        $sum{$k} = $#{$map{$k}};
        for my $i (0..$#{$map{$k}}){
            if ($_ eq @{$map{$k}}[$i]) {
                $count{$k}++;
                @{$map{$k}}[$i] .= "(m)";
            }
        }
    }
}

for my $i (keys %map) {
    print "$i\t$count{$i}\/$sum{$i}\t@{$map{$i}}\n";
}

