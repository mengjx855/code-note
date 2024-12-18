#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-03-05, 10:38:33
# modified date: 2024-02-10, 21:20:09
shopt -s expand_aliases

if [ $# -lt 4 ];then
    echo "$0 [fq|fq1,fq2] [contig] [prefix] [out_directory] [min_contig:-1500]"
    exit 127
fi
# 输入out_directory，输出bin文件位于该文件夹下bin_fa目录，中间文件位于该文件夹下tmp_binning目录

fq=$1
fa=$2
p=$3
tmp=$4/tmp_binning
out=$4/bin_fa
len=${5:-1500}
trds=12

alias bwa-mem2='/share/data1/software/bwa-mem2-2.2.1/bwa-mem2'
alias sam2dep='/share/data1/mjx/scripts/tools_get_depth_from_sam.single.py'
alias metabat2='/share/data1/software/miniconda3/envs/metabat2/bin/metabat2'

( [ ! -d $tmp ] ) && mkdir -p $tmp
( [ -e $tmp/$p.log ] ) && { echo -e "Skip sample: $p .." && exit 0; }

if [ ! -e $tmp/$p.depth ];then 
    bwa-mem2 index -p $tmp/$p.index $fa >> $tmp/$p.log 2>&1 || { rm $tmp/$p.log $tmp/$p.index* && exit 1; }
    if [[ $fq =~ "," ]];then
        fq1=$(echo $fq | cut -d "," -f1)
        fq2=$(echo $fq | cut -d "," -f2)
        bwa-mem2 mem $tmp/$p.index $fq1 $fq2 -t $trds 2>> $tmp/$p.log | sam2dep $p > $tmp/$p.depth || { rm $tmp/$p.depth $tmp/$p.log && exit 1; }
    else
        bwa-mem2 mem $tmp/$p.index $fq -t $trds 2>> $tmp/$p.log | sam2dep $p > $tmp/$p.depth || { rm $tmp/$p.depth $tmp/$p.log && exit 1; }
    fi
    rm $tmp/$p.index*
    chmod 444 $tmp/$p.depth
fi

( [ ! -d $out ] ) && mkdir -p $out
metabat2 -t $trds -m $len -s 200000 --saveCls --unbinned --seed 2024 -i $fa -a $tmp/$p.depth -o $out/$p.bin >> $tmp/$p.log ||\
    { rm $tmp/$p.log $out/$p.bin* && exit 1; }
chmod 444 $tmp/$p.log $out/$p.bin*
mv $out/$p.bin $tmp/$p.bin.Cls >/dev/null 2>/dev/null
mv $out/$p.bin.unbinned.fa $out/$p.bin.lowDepth.fa $out/$p.bin.tooShort.fa $tmp >/dev/null 2>/dev/null
