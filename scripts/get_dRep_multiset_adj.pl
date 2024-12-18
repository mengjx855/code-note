#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl $0 <dRep_cls2_info> <sig_list> <out_f>\nInfo: sig must be unique in certain genome set, and should supply character in a comma-delimited format.\n" unless @ARGV eq 3;
my ($f, $sig ,$o) = @ARGV;

open IN, "<$f" or die "Can't open $f: $!\n";
open OU, ">$o"; 

my @s = split /,/, $sig;
print OU "cls\tsig_".join("\tsig_", @s)."\n";

while(<IN>){
    chomp;
    my @s2 = split /\t/, $_;
    my @s3 = split /;/, $s2[2];
    
    my %h;
    foreach my $i (@s) {$h{$i} = 0;}
    foreach my $i (@s3) {
        foreach my $pattern (@s) {
            if ($i =~ /$pattern/) {
                $h{$pattern} = $h{$pattern} + 1;
            } else { next }
        }
    }
    my @l;
    foreach my $i (@s) {
        push @l, $h{$i};
    }
    print OU "$s2[0]\t".join("\t", @l)."\n";
}
print STDERR "Program End...\n";

###########################################################
# - Meng Jinxin
# - 2023-05-08