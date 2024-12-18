#!/usr/bin/bash
# Author: Jinxin Meng
# Created Date: 2023-01-01
# Modified Date: 2023-09-30
# Version: 1.0

if [ $# -ne 2 ];then
    echo "$0 <rc> <rc_sum>"
    exit 127
fi

rc=$1
rc_sum=$2

perl -e '$sum=0;while(<>){chomp;@l=split/\s+/;$l[0]=~/(.*)_\d+$/;if ($pre eq $1){$sum+=$l[1]}else{if ($pre ne ""){print "$pre\t$sum\n"};$sum=$l[1];$pre=$1}};print "$pre\t$sum\n"' $rc > $rc_sum
