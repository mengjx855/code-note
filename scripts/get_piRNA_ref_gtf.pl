#!/share/software/anaconda3/bin/perl
use warnings;
use strict;

################################################
### Author: mjx/mengjx855@163.com            ###
### Created Date: 2022-07-09                 ###
### Version: 1.0                             ###
################################################

die "Usage: perl $0 <genome.fasta> <piR_gtf> <piR_ref_gtf> \n" unless @ARGV == 3;
my ($genome, $piR_gtf, $piR_ref_gtf) = @ARGV;

my %h = ();

###get genomic chr id###
open IN, "<$genome" or die $!;
while(<IN>){
    chomp;
    if(/>/){
        my @a = split /\s+/;
        $a[0] =~ s/>//;
        $h{$a[1]} = $a[0];
    }
}
close IN;
###

###piR_gtf2piR_ref_gtf###
open IN, "<$piR_gtf" or die $!;
open OU, ">$piR_ref_gtf" or die $!;
while(<IN>){
    chomp;
    if(/^\#/){
       print OU "$_\n";
   } else {
       my @s = split /\t/;
       $s[2] = "exon";
       $s[8] =~ m/piRNA_id \"piR-hsa-(\d+)\";/;
       my $piRid = $1;
       my $old = shift @s;
       delete $s[7];
       my $piRinfo = "gene_id \"PIRG$piRid\"; transcript_id \"PIRT$piRid\"; gene_type \"piRNA\"; gene_status \"KNOWN\"; gene_name \"piR-hsa-$piRid\"; transcript_type \"piRNA\"; transcript_status \"KNOWN\"; transcript_name \"piR-hsa-$piRid\"; exon_id \"PIRE$piRid\"; ";
       print OU "$h{$old}\t".join("\t",@s)."\t$piRinfo\n";
   }
}
close IN;
close OU;
###
print STDERR "Program End...\n";
