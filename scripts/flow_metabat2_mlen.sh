#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-03-05, 10:38:33
# modified date: 2024-02-10, 21:20:09
shopt -s expand_aliases

if [ $# -lt 4 ];then
    echo "$0 <fq|fq1,fq2> <contig> <prefix> <out_directory> <min_contig:-1500>"
    exit 127
fi
# 输入out_directory，输出bin文件位于该文件夹下bin_fna目录，中间文件位于该文件夹下tmp_binning目录

fq=$1
fa=$2
p=$3
tmp=$4/tmp_binning
out=$4/bin_fna
len=${5:-1500}
trds=24

alias seqkit='/share/data1/software/binary/seqkit'
alias bwa='/share/data1/software/bwa-0.7.17/bwa'
alias sam2dep='/share/data1/mjx/bin/tools_get_depth_from_sam.single.py'
alias jgi_summarize_bam_contig.depths='/share/data1/software/miniconda3/envs/metabat2/bin/jgi_summarize_bam_contig.depths'
alias metabat2='/share/data1/software/miniconda3/envs/metabat2/bin/metabat2'

if [ ! -d $tmp ];then 
    mkdir -p $tmp
fi

if [ -e $tmp/$p.log ];then
    echo -e "Skip sample: $p .."
    exit 0
fi

if [ ! -e $tmp/$p.m$len.fna ];then
    seqkit seq -g --quiet -m $len $fa -o $tmp/$p.m$len.fna || { rm $tmp/$p.m$len.fna && exit 1; }
fi

if [ ! -e $tmp/$p.depth ];then 
    bwa index -p $tmp/$p.index $tmp/$p.m$len.fna 2>> $tmp/$p.log || { rm $tmp/$p.log $tmp/$p.index* && exit 1; }
    if [[ $fq =~ "," ]];then
        fq1=$(echo $fq | cut -d "," -f1)
        fq2=$(echo $fq | cut -d "," -f2)
        bwa mem $tmp/$p.index $fq1 $fq2 -t $trds 2>> $tmp/$p.log | sam2dep $p > $tmp/$p.depth && chmod 444 $tmp/$p.depth || { rm $tmp/$p.log && exit 2; }
    else
        bwa mem $tmp/$p.index $fq -t $trds 2>> $tmp/$p.log | sam2dep $p > $tmp/$p.depth && chmod 444 $tmp/$p.depth || { rm $tmp/$p.log && exit 2; }
    fi
    rm $tmp/$p.index*
fi

if [ ! -d $out ];then
    mkdir -p $out
fi

metabat2 -t $trds -m $len -s 200000 --saveCls --unbinned --seed 2024 -i $tmp/$p.m$len.fna -a $tmp/$p.depth -o $out/$p.bin >> $tmp/$p.log &&\
    chmod 444 $tmp/$p.log || { rm $tmp/$p.log && exit 2; }
mv $out/$p.bin $tmp/$p.bin.Cls >/dev/null 2>/dev/null
mv $out/$p.bin.unbinned.fa $out/$p.bin.lowDepth.fa $out/$p.bin.tooShort.fa $tmp >/dev/null 2>/dev/null
rm $tmp/$p.m$len.fna