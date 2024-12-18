#!/usr/bin/bash
# Author: Jinxin Meng
# Created Date: 2023-07-27
# Modified Date: 2023-09-19
# Version: 1.0

if [ $# -ne 2 ];then 
    echo -e "Usage: bash $0 <configs> <out_directory>
\e[1;31m Configs format (tab-delimited file): \e[0m
    <prefix>    <fastq1>    <fastq2>    <contigs>
    s1|..       s1_1.fq.gz  s1_2.fq.gz  s1.contigs.fa
    s2|..       s2_1.fq.gz  s2_2.fq.gz  s2.contigs.fa
    s3|..       s3_1.fq.gz  s3_2.fq.gz  s3.contigs.fa
    s4|..       s4_1.fq.gz  s4_2.fq.gz  s4.contigs.fa
 Binning (metabat2) multiple coverage mode was performed between s1, s2, s3, s4."
    exit 127
fi

cfg=$1
out=$2
tmp=$(dirname $out)/tmp_binning
len=2000
jobs=5
trds=16

alias parallel=/usr/local/bin/parallel
alias seqkit=/share/data1/software/binary/seqkit
alias bwa=/share/data1/software/bwa-0.7.17/bwa
alias samtools=/share/data1/software/samtools/bin/samtools
alias jgi_summarize_bam_contig_depths=/share/data1/software/miniconda3/envs/metabat2/bin/jgi_summarize_bam_contig_depths
alias metabat2=/share/data1/software/miniconda3/envs/metabat2/bin/metabat2
alias bzip2=/share/data1/software/miniconda3/bin/bzip2
alias binning_dep_combine.pl=/share/data1/mjx/bin/binning_dep_combine.pl

# print_log(){ 
#     echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] $*" 
# }

declare -A ctg
declare -A fq
while read i;do
    line=($(echo $i | awk -F '\t' '{print $1 $2 $3 $4}'))
    ctg[${line[0]}]=${line[3]}
    fq[${line[0]}]=`echo -e "${line[1]}\t${line[2]}"`
done < $cfg
# echo "Array key: ${!ctg[*]}"

if [ ! -d $tmp ];then mkdir -p $tmp;fi
for i in ${!ctg[*]};do
    if [ ! -e ${tmp}/${i}.contigs.m${len}.fa ];then echo -e "$i\t${ctg[$i]}";fi
done | parallel -j ${jobs} --colsep="\t" /usr/bin/time -v -a -o ${tmp}/time.log seqkit seq -g -m ${len} -o ${tmp}/{1}.contigs.m${len}.fa {2}
# print_log "seqkit: contigs length less than ${len} were filtered completion.." ${tmp}/metabat2_mcvg.log 

for i in ${!ctg[*]};do
    if [ ! -e ${tmp}/${i}.contigs.m${len}.bwt ];then echo -e "${i}";fi
done | parallel -j ${jobs} /usr/bin/time -v -a -o ${tmp}/time.log bwa index -p ${tmp}/{} ${tmp}/{}.contigs.m${len}.fa 2\>\> ${tmp}/{}.log
# print_log "bwa build completion.." ${tmp}/metabat2_mcvg.log

for i in ${!ctg[*]};do
    if [ ! -e ${tmp}/${i}.depth ];then
        for j in ${!fq[*]};do
            echo -e "$i\t$j\t${fq[${j}]}"
        done | \
            parallel -j ${jobs} --colsep="\t" /usr/bin/time -v -a -o ${tmp}/time.log \
                bwa mem ${tmp}/{1} {3} {4} -t ${trds} -o ${tmp}/{2}_map_to_{1}.sam 2\>\> ${tmp}/{1}.log \&\& \
                /usr/bin/time -v -a -o ${tmp}/time.log samtools view -@ ${trds} -bS ${tmp}/{2}_map_to_{1}.sam \| \
                    samtools sort -@ ${trds} -o ${tmp}/{2}_map_to_{1}_sort.bam \> /dev/null 2\>1 \&\& \
                rm ${tmp}/{2}_map_to_{1}.sam

        rm ${tmp}/${i}.amb ${tmp}/${i}.ann ${tmp}/${i}.bwt ${tmp}/${i}.pac ${tmp}/${i}.sa
        
        for j in ${!fq[*]};do
            echo -e "$i\t$j"
        done | \
            parallel -j ${jobs} --colsep="\t" /usr/bin/time -v -a -o ${tmp}/time.log \
                jgi_summarize_bam_contig_depths --outputDepth ${tmp}/{2}_map_to_{1}.depth ${tmp}/{2}_map_to_{1}_sort.bam 2\>\> ${tmp}/{1}.log \&\& \
                rm ${tmp}/{2}_map_to_{1}_sort.bam
        binning_dep_combine.pl ${tmp}/${i}.contigs.m${len}.fa ${tmp}/*_map_to_${i}.depth > ${tmp}/${i}.depth
        bzip2 ${tmp}/*_map_to_${i}.depth
    fi
done
# print_log "bwa mem, samtools view and sort, jgi_summarize_bam_contig_depths calu depth, completion.." >> ${tmp}/metabat2_mcvg.log

if [ ! -d ${out}_fna ];then mkdir ${out}_fna;fi
if [ ! -d ${out}_unbinned ];then mkdir ${out}_unbinned;fi
for i in ${!ctg[*]};do
    echo -e "$i\t${tmp}/${i}.contigs.m${len}.fa\t${tmp}/${i}.depth"
done | \
    parallel -j ${jobs} --colsep="\t" /usr/bin/time -v -a -o ${tmp}/time.log \
        metabat2 -t ${trds} -m ${len} -s 200000 --saveCls --unbinned --seed 2023 \
            -i {2} -a {3} -o ${out}_fna/{1}.bin \>\> ${tmp}/{1}.log \&\& \
        bzip2 ${tmp}/{1}.depth \&\& \
        mv ${out}_fna/{1}.bin ${tmp}/{1}.bin.cls \&\& \
        mv ${out}_fna/{1}.bin.unbinned.fa ${out}_fna/{1}.bin.lowDepth.fa ${out}_fna/{1}.bin.tooShort.fa ${tmp} \&\& \
        rm ${tmp}/{1}.contigs.m${len}.fa
# print_log "metabat2 binning completion.." >> metabat2_mcvg.log

chmod 444 ${tmp}/*.log ${tmp}/*.depth*
# print_log "Program End.." >> metabat2_mcvg.log

