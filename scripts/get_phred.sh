#!/usr/bin/bash
# encoding: utf-8
# author : Jinxin Meng
# created date: 2022-01-02, 15:21:49
# modified date: 2023-12-02, 15:21:49
if [ $# -lt 1 ];then
    echo "$0 <fq|fq.gz>"
    exit 127
fi

fq=$1

less $fq |\
    head -n 10000 |\
    awk '{if(NR%4==0) printf("%s",$0);}' |\
    od -A n -t u1 -v |\
    awk 'BEGIN{min=100;max=0;} \
    {for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END \
        {if(max<=126 && min<59) print "'$fq' Phred33"; \
        else if(max>73 && min>=64) print "'$fq' Phred64"; \
        else if(min>=59 && min<64 && max>73) print "'$fq' Solexa64"; \
        else print "'$fq' Unknown score encoding"}'