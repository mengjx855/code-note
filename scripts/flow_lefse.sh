#!/usr/bin/bash

set -e # 如果出错，就不再执行下一步

if [ $# -lt 2 ];then
    echo "$0 <input_otu> <output_prefix>"
    exit 127
fi

in_f=$1
out_f=$2

source activate /home/mengjx/.conda/envs/lefse

lefse_format_input.py $in_f  ${out_f}.fmt -c 1 -o 1000000

lefse_run.py ${out_f}.fmt ${out_f}.res  -l 2
#awk '$4!="-" && $3!="-"' ${out_f}.res | perl -ne 'chomp;@l=split/\s+/;@x=split(/\./, $l[0]);$n=join(".", @x[])'

lefse_plot_cladogram.py ${out_f}.res ${out_f}.cladogram.pdf --format pdf \
            --max_lev 7 --labeled_start_lev 2 --labeled_stop_lev 7 --abrv_start_lev 3 --abrv_stop_lev 7 \
                        --right_space_prop 0.3 --title "log10 (LDA) >2"

lefse_plot_res.py --format pdf ${out_f}.res ${out_f}.pdf
