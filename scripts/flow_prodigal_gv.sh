#!/usr/bin/bash

set -e

if [ $# -lt 2 ] || [ $# -gt 4 ];then
echo "$0 <fasta> <out_prefix> <threads> <translate_tab|11>"
echo "        <out_prefix>.faa/gff/ffn"
exit 127
fi
in_f=$1
out_f=$2
threads=$3
translate_tab=${4:-11}
temp_prodital=${out_f}_$$  # $$ is pid of this script 

mkdir $temp_prodital

parallel -k  -j ${threads} -a ${in_f} --block -1 \
             --pipe-part --recend '\n' --recstart '>' \
             "cat | \
             prodigal-gv -q -d $temp_prodital/{#}.ffn -o $temp_prodital/{#}.gff  -p meta -a $temp_prodital/{#}.faa -f gff -g ${translate_tab}"
cat $temp_prodital/*.faa > $out_f.faa
cat $temp_prodital/*.ffn > $out_f.ffn
cat $temp_prodital/*.gff > $out_f.gff
rm -r $temp_prodital
