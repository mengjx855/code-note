#!/usr/bin/bash
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-06-07, 17:02:49
# modified date: 2024-01-22, 10:15:55

if [ $# -ne 2 ];then
    echo "$0 <gene_faa> <gene_gff> <out directory>"
    exit 127
fi

fas=$1
gff=$2
out=$3
trds=80
db=/share/data1/database/run_dbcan4/db_20231207

source activate /share/data1/software/miniconda3/envs/run_dbcan
run_dbcan --db_dir $db --out_dir $out \
    --dia_cpu $trds --dia_eval 1e-102 \
    --hmm_cpu $trds --hmm_cov 0.35 --hmm_eval 1e-15 \
    --dbcan_thread $trds \
    --tf_cpu $trds --tf_eval 0.00001 --tf_cov 0.35 \
    --stp_cpu $trds --stp_eval 0.00001 --stp_cov 0.3 \
    --cluster $gff --cgc_dis 2 --cgc_sig_genes all \
    --cgc_substrate \
    $fas protein
conda deactivate
