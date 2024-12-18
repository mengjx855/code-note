#!/usr/bin/perl
# Author: Jinxin Meng
# Create Date: 2023-80-04
# Modified Date: 2023-10-07
# Source: WatsonLab/single_and_multiple_binning/scripts/combine.pl
# Version: 1.0

use strict;
use File::Basename;

die "Usage: perl $0 [contigs] jgiSummarizeBamContigDepths1 [ jgiSummarizeBamContigDepths2 ...]\n" unless @ARGV gt 1;

my $ctg = $ARGV[0];
shift @ARGV;
my @files = @ARGV;
my $o;
my @names;
my $tadcounter = 0;

foreach my $file (@files) {
    open(IN, $file);
    while(<IN>) { last }

    $tadcounter++;
    my $base = basename $file, ".depth";
    $base = $base . ".bam";
    push(@names, $base);

    while(<IN>) {
        chomp();
        my($cn, $cl, $tad, $bam, $var) = split(/\t/);
        $o->{$cn}->{contigInfo}->{cl} = $cl;
        $o->{$cn}->{contigInfo}->{tad} += $tad;
        $o->{$cn}->{bamInfo}->{$base}->{bam} = $bam;
        $o->{$cn}->{bamInfo}->{$base}->{var} = $var;
    }
    close IN;
}

print "contigName\tcontigLen\ttotalAvgDepth";
foreach my $name (@names) {
    print "\t$name\t$name" . "-var";
}
print "\n";

open IN, "<$ctg";
my @cnlist;
while(<IN>){
    chomp;
    if(/>(\S+)/){
        push @cnlist, $1;
    }
}
close IN;

foreach my $cn (@cnlist) {

    my $cr = $o->{$cn}->{contigInfo};
    my $br = $o->{$cn}->{bamInfo};
    print $cn, "\t", $cr->{cl}, "\t", $cr->{tad} / $tadcounter;

    foreach my $name (@names) {
        print "\t", $br->{$name}->{bam}, "\t", $br->{$name}->{var};
    }
    print "\n";
}

=pod
sub sortfunction {
    my($a1,$a2) = split(/_/, $a);
    my($b1,$b2) = split(/_/, $b);
    $a1 =~ s/^k//;
    $b1 =~ s/^k//;
    return $a1 <=> $b1 || $a2 <=> $b2;
}
=cut
