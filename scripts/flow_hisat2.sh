#!/usr/bin/bash
# encoding: utf-8
# author : Jinxin Meng
# created date: 2022-03-21, 22:54:38
# modified date: 2023-11-18, 22:54:38

if [ $# -lt 3 ];then
    echo "$0 <fq|fq1,fq2> <index> <out_prefix|sort.bam>"
    echo -e "    \e[1;31mFor example:\e[0m\n For single-end: ${0##*/} prefix.fq.gz index prefix\n For pair-end: ${0##*/} prefix_1.fq.gz,prefix_2.fq.gz index prefix"
    echo -e "    \e[1;31mCurrent available hisat2 index:\e[0m
 None"
    exit 127
fi

fq=$1
idx=$2
p=$3
trds=16
alias hisat2=/share/software/hisat2-2.2.1/hisat2
alias samtools=/share/software/binary/samtools

if [ -f ${p}_hisat2.log ] && grep -q 'overall' ${p}_hisat2.log;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] INFO Existing and skip sample: $p."
    exit 0
fi

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    hisat2 -p $trds -x $idx -1 $fq1 -2 $fq2 2>${p}_hisat2.log |\
        samtools view -@ $trds -bS - |\
        samtools sort -@ $trds -o ${p}_sort.bam - \>/dev/null 2\>/dev/null
else
    hisat2 -p $trds -x $idx -U $fq 2>${p}_hisat2.log |\
        samtools view -@ $trds -bS - |\
        samtools sort -@ $trds -o ${p}_sort.bam - \>/dev/null 2\>/dev/null
fi
chmod 444 ${p}_hisat2.log ${p}_sort.bam
