#!/usr/bin/bash

set -e # 如果出错，就不再执行下一步

if [ $# -lt 2 ];then
    echo "$0 <input> <output>"
    exit 127
fi

input=$1
output=$2

source activate /home/mengjx/.conda/envs/mengjx/
mkdir -p $output/00.cor/ $output/01.bsp/ $output/02.bsp_cor/ $output/03.pval/

python /share/data1/mengjx/software/git/SparCC3/SparCC.py $input -i 4 -c $output/00.cor/sparcc_cor.txt -v $output/00.cor/sparcc_cov.txt >> $output/00.cor/sparcc.log
python /share/data1/mengjx/software/git/SparCC3/MakeBootstraps.py $input -n 100 -t bootstrap_#.txt -p $output/01.bsp/ >> $output/01.bsp/sparcc.log
# for n in {0..99}; do python /share/data1/mengjx/software/git/SparCC3/SparCC.py $output/01.bsp/bootstrap_${n}.txt -i 10 -c $output/02.bsp_cor/bootstrap_cor_${n}.txt -v $output/02.bsp_cor/bootstrap_cov_${n}.txt >> $output/02.bsp_cor/sparcc.log; done
for n in {0..99};do echo python /share/data1/mengjx/software/git/SparCC3/SparCC.py $output/01.bsp/bootstrap_${n}.txt -i 4 -c $output/02.bsp_cor/bootstrap_cor_${n}.txt -v $output/02.bsp_cor/bootstrap_cov_${n}.txt \>\> $output/02.bsp_cor/sparcc.log;done | parallel -j 20 {}
# seq 0 1 99 | parallel -j 30 python /share/data1/limh/software/SparCC3-master/SparCC.py $output/01.bsp/bootstrap_{}.txt -i 4 -c $output/02.bsp_cor/bootstrap_cor_{}.txt -v $output/02.bsp_cor/bootstrap_cov_{}.txt >> $output/02.bsp_cor/sparcc.log &
python /share/data1/mengjx/software/git/SparCC3/PseudoPvals.py $output/00.cor/sparcc_cor.txt $output/02.bsp_cor/bootstrap_cor_#.txt 100 -o $output/03.pval/sparcc_pvals.txt -t two_sided >> $output/03.pval/sparcc.log

cd $output;ln -s 00.cor/sparcc_cor.txt;ln -s 03.pval/sparcc_pvals.txt
conda deactivate
