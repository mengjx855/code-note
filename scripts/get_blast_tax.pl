#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-02-22, 12:18:19
# modified date: 2024-03-08, 10:29:48
use warnings;
use strict;
use Getopt::Long;

my ($blast, $db, $out_f, $fmt, $help);
GetOptions(
    "i|in:s" => \$blast,
    "d|db:s" => \$db,
    "o|out:s" => \$out_f,
    "out_blast+" => \$fmt,
    "h|help:s" => \$help,
);

&usage if (! $blast or ! $db or ! $out_f or $help);

my (%h, %ref);
$h{"prok"} = '/share/data1/database/ncbi/nt_prok_20230924/nt_prok_accession2tax';
$h{"virus"} = '/share/data1/database/ncbi/nt_viruses_20230924/nt_viruses_accession2tax';
$h{"fungi"} = '/share/data1/database/ncbi/nt_euk_fungi_20240307/nt_euk_fungi_accession2tax';
open I, "$h{$db}" or die "Can't open $h{$db}: $!\n";
while (<I>) {
    chomp;
    my @s = split/\s+/;
    next if $#s eq 1;
    $ref{$s[0]} = "$s[1]\t$s[2]";
}

open I, "$blast" or die "Can't open $blast: $!\n";
open O, ">$out_f" or die "Can't open $out_f: $!\n";
while (<I>) {
    chomp;
    my @s = split/\s+/;
    if (exists $ref{$s[1]}) {
        if ($fmt) {
            print O "$_\t$ref{$s[1]}\n";
        } else {
            print O "$s[0]\t$ref{$s[1]}\n";
        }
    } else {
        if ($fmt) {
            print O "$_\tNULL\tNULL\n";
        } else {
            print O "$s[0]\tNULL\tNULL\n";
        }
    }
}
close O;

sub usage {
    print  "\e[;32;1m\n","Usage: perl $0 -i blast_out -d *accession2tax -o out_file\n","\e[;32;1m\n";
    die "Options
    This program is used to obtain taxonomic info. for blast out.
    -i/--in FILE    input blast out with format 6
    -d/--db FILE    reference accession2tax with field \"asscession|taxid|info.\" [prok, viruses, fungi]
    -o/--out FILE   output tab-delimited file
    --out_blast     blast-like output
    -h/--help       print usage and DESCRIPTION, ignore all other parameters
    \n";
}

