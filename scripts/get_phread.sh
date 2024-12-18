#!/usr/bin/bash
# encoding: utf-8
# author : Jinxin Meng
# created date: 2022-01-02, 15:21:49
# modified date: 2024-02-05, 18:08:44
if [ $# -lt 1 ];then
    echo "$0 <fq|fq.gz>"
    exit 127
fi
less $1 |\
    head -n 10000 |\
    awk '{if(NR%4==0) printf("%s",$0);}' |\
    od -A n -t u1 -v |\
    awk 'BEGIN{min=100;max=0;} \
    {for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END \
        {if(max<=126 && min<59) print "'$1' Phred33"; \
        else if(max>73 && min>=64) print "'$1' Phred64"; \
        else if(min>=59 && min<64 && max>73) print "'$1' Solexa64"; \
        else print "'$1' Unknown score encoding"}'
