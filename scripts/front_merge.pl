#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-05-31, 09:34:08
# modified date: 2024-02-19, 09:34:08
use warnings;
use strict;

die "perl $0 [in_f1] [in_f2] [..] [out_file]\n" unless @ARGV ge 3;
my ($in_f, $out_f) = ($ARGV[0], $ARGV[-1]);
my (@list, %h, %h2) = ();

open O,">$out_f";
pop @ARGV;
my $len = @ARGV;
foreach (@ARGV) {
    open I,"<$_";
    while(<I>){
        chomp;
        my @s = split/\t/;   
        $h{$s[0]}++;
    }
    close I;
}

foreach (keys %h) { 
    $h2{$_}++ if $h{$_} eq $len;
}

open I,"<$in_f";
while(<I>) {
    chomp;
    if (/contig_id/) {
        print O "$_\n";
        next
    }
    my @s = split/\t/;
    print O "$_\n" if exists $h2{$s[0]}; 
}
close O;

