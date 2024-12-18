#!/usr/bin/bash
# Author: Jinxin Meng
# Created Date: 2023-01-01
# Modified Date: 2023-09-30
# Version: 1.0

if [ $# -lt 3 ];then
    echo "$0 <fq|fq1,fq2> <index> <out_prefix> <u:-10000000>"
    exit 127
fi

fq=$1
index=$2
p=$3
u=${4:-10000000}
trds=16

alias bowtie2='/share/data1/software/bowtie2-2.5.1/bowtie2'

if [[ ${fq} =~ "," ]] ; then
    fq1=$(echo ${fq} | cut -d "," -f1)
    fq2=$(echo ${fq} | cut -d "," -f2)

    bowtie2 --end-to-end --mm --no-head --no-unal --no-sq \
        -u $u -1 ${fq1} -2 ${fq2} -x ${index} -S ${p}.sam -p ${trds} 2> ${p}.log
else
    bowtie2 --end-to-end --mm --no-head --no-unal --no-sq \
        -u $u -U ${fq} -x $index -S ${p}.sam -p ${trds} 2> ${p}.log
fi
perl -e 'while(<>){if(/^\S+\s+\S+\s+(\S+)\s+/){$h{$1}++;}} for(sort keys %h){print "$_\t$h{$_}\n";}' ${p}.sam > ${p}.rc
rm ${p}.sam
