#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-12-13, 11:37:12
# modified date: 2023-12-13, 13:21:21

if [ $# -lt 2 ]
then
    echo "$0 <contig_fna> <out_prefix> <out_directory> <ani:95> <cvg:70> <len_len:3000> <threads:32>"
    echo "  Readme: parameter of ani, cvg, len_len and threads could be modified in the script if you need."
    exit 2
fi

ctg=$1
p=$2
out=${3:-.}
ani=95
cvg=70
len=2000
trds=32
alias makeblastdb='/share/data1/software/ncbi-blast-2.13.0+/bin/makeblastdb'
alias blastn='/share/data1/software/ncbi-blast-2.13.0+/bin/blastn'
alias anicalc.py='/share/data1/software/ckv_scripts/anicalc.py'
alias aniclust.py='/share/data1/software/ckv_scripts/aniclust.py'
alias seqkit='/share/data1/software/binary/seqkit'

if [ ! -d $out ];then
    mkdir $out
fi

if [ -e $out/$p.log ];then
    echo -e "Skip sample: $p .."
    exit 0
fi
seqkit seq -g --quiet -m $len $ctg -o $out/${p}_m${len}.fna || { rm $out/${p}_m${len}.fna && exit 1; }
makeblastdb -in $out/${p}_m${len}.fna -out $out/${p}_m${len}_idx -dbtype nucl -hash_index >> $out/$p.log ||
    { rm $out/${p}_m${len}_idx* $out/$p.log && exit 1; }
blastn -query $out/${p}_m${len}.fna -db $out/${p}_m${len}_idx -out $out/${p}_m${len}_btn -outfmt '6 std qlen slen' \
    -num_alignments 99999999 -num_threads $trds >> $out/$p.log || { rm $out/${p}_m${len}_btn $out/$p.log && exit 1; }
anicalc.py -i $out/${p}_m${len}_btn -o $out/${p}_m${len}_ani.tsv >> $out/$p.log || { rm $out/${p}_m${len}_ani.tsv $out/$p.log && exit 1; }
aniclust.py --fna $out/${p}_m${len}.fna --ani $out/${p}_m${len}_ani.tsv --out $out/${p}_m${len}_cluster.tsv \
    --min_ani $ani --min_tcov $cvg --min_qcov 0 >> $out/$p.log || { rm $out/${p}_m${len}_cluster.tsv $out/$p.log && exit 1; }
cut -f1 $out/${p}_m${len}_cluster.tsv | seqkit grep --quiet -f - -o $out/${p}_m${len}_cluster.fna $out/${p}_m${len}.fna ||\
    { rm $out/${p}_m${len}_cluster.fna $out/$p.log && exit 1; }
rm $out/${p}_m${len}_idx* $out/${p}_m${len}.fna
chmod 444 $out/${p}_m${len}_btn $out/${p}_m${len}_cluster.fna $out/$p.log $out/${p}_m${len}_ani.tsv $out/${p}_m${len}_cluster.tsv
