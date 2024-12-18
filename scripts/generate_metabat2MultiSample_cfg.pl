#!/usr/bin/perl
use strict;

die usage()."\n" unless @ARGV eq 3; 

my ($sampleSet, $ctgDir, $fqDir) = @ARGV;
my @samples;

$ctgDir = $ctgDir."\/" unless $ctgDir =~/\/$/;
$fqDir = $fqDir."\/" unless $fqDir =~/\/$/;

open IN, "<$sampleSet";
while(<IN>){chomp; push @samples, $_;}

foreach my $sample (@samples) {
    # print sampleID, contigDir, fqDir
    print "$sample\t$ctgDir$sample.contigs.fa\t$fqDir$sample\_dehost_1.fq.gz\t$fqDir$sample\_dehost_2.fq.gz\n";
}

sub usage {
    return "Usage: perl $0 [SampleSet] [contigDirectory] [Paired-End fastqDirectory]\nDefault:
    contigs suffix: [sample]_contigs.fa;
    fastq suffix: [sample]_dehost_[1/2].fq.gz."
}

print STDERR "Program End...\n";
###########################################################
## Source: generate_metabat2MultiSample_cfg.pl
## Author: Jinxin Meng
##   Data: 2023-07-25