#!/usr/bin/bash

set -e # 如果出错，就不再执行下一步

if [ $# -lt 4 ];then
   echo "$0 <fq1> <fq2> <ref> <prefix>"
   exit 127
fi

fq1=$1
fq2=$2
ref=$3
prefix=$4

bowtie2-build $ref $ref

bowtie2 --end-to-end --sensitive --no-unal -p 16 -x $ref -1 $fq1 -2 $fq2 --un-conc-gz $prefix -S $prefix.sam
