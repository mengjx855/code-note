#!/usr/bin/bash
##########################################################
# Creater       :  夜下凝月
# Created  date :  2023-09-08, 16:42:31
# Modiffed date :  2023-09-08, 16:42:31
##########################################################


( [ $# -ne 3 ] ) && echo "$0 <fasta> <depth> <path_output_prefix>" && exit 127

in_fa=`realpath $1` # file
depth=`realpath $2` # file
pref=`realpath $3` # prefix of output. (relative/absolute path)

pref_name=${pref##*/} # 输出文件的名字

## 如果已经有了bins，那就退出
( [ -f ${pref}.bins.ok ] ) && echo "exists ${pref}.bins" && exit 127

metabinner_dir="/share/data1/software/MetaBinner-master"

( [ ! -d ${pref} ] ) && { mkdir -p "${pref}" || exit 127 ;}
cd ${pref}

source activate /home/zhangy/.conda/envs/metabinner_env  || exit 127

## link fasta
( [ ! -L ${pref_name}.tmp.fa ] || [ ! -f ${pref_name}.tmp.fa ] ) && { ln -s  ${in_fa} ${pref_name}.tmp.fa || exit 127 ;}
## format depth
( [ ! -f ${pref_name}.tmp.depth ] ) && { cut -f 1,4 ${depth} > ${pref_name}.tmp.depth || exit 127 ;}
## generate kmer
( [ ! -f ${pref_name}.kmer.ok ] ) && { python ${metabinner_dir}/scripts/gen_kmer.py.zy ${pref_name}.tmp.fa 0 4 ${pref_name}.kmer.csv && touch ${pref_name}.kmer.ok || ! rm ${pref_name}.kmer.csv || exit 127 ;}
## get cluster
( [ ! -f ${pref_name}.logtrans.tsv ] ) && { python ${metabinner_dir}/scripts/component_binning.py.zy \
    --software_path ${metabinner_dir} \
    --contig_file ${pref_name}.tmp.fa \
    --coverage_profiles ${pref_name}.tmp.depth \
    --composition_profiles ${pref_name}.kmer.csv \
    --output ./ --log ${pref_name}.out.log \
    --threads 80 \
    --dataset_scale large \
    --contig_length_threshold 500 \
    && touch ${pref_name}.logtrans.ok \
    || exit 127 ;}

## get bins
( [ ! -f ${pref_name}.bin.ok ] ) && { 
python ${metabinner_dir}/scripts/gen_bins_from_tsv.py \
    -f ${pref_name}.tmp.fa \
    -r ${pref}/intermediate_result/kmeans_length_weight_X_t_logtrans_result.tsv \
    -o ${pref_name}.bins \
    && mv ${pref_name}.bins ${pref}.bins \
    && touch ${pref_name}.bin.ok \
    && cd .. \
    && touch ${pref}.bins.ok \
    || exit 127 ;}

## remove tmp
( [ -f ${pref_name}.bins.ok ] ) && { rm -rf ${pref} || exit 127; }
