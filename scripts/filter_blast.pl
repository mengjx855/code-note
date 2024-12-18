#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-06-06, 11:42:16
# modified date: 2024-06-06, 11:42:16

use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage .= "Usage: perl $0 [OPTIONS] [*.btn or *.btp] [out_file]\n";
$usage .= "  OPTIONS:\n"
$usage .= "   --top1 : select first items\n";
$usage .= " --by_bit : sort by bitscore\n";
$usage .= " --by_cvg : sort by coverage\n";
$usage .= " --by_idt : sort by identity\n";


my ($ckm, $ckm2, $len) = ("", "", "");

&GetOptions(
    "ckm=s" => \$ckm,
    "ckm2=s" => \$ckm2,
    "len=s" => \$len);

die "$usage" unless @ARGV eq 2;
my ($cdb, $out) = @ARGV;

