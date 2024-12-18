#!/usr/bin/bash

if [ "$#" -lt "7" ]; then
  echo -e "\nusage: sh $0 <Map|SNP|All> <fastq F> <fastq R> <database prefix(for bowtie2 and mpileup)> <database ffn> <output file name> <output dir> \n" 
    exit 2
fi

Type=$1
fq1=$2
fq2=$3
db=$4
gffPath=$5
Samp=$6
outdir=$7

readPair=1000000000
threads=13
SNP_cov=0.05
SNP_depth=40
min_site=200000
lab='[.]'

#相关应用程序路径
#export PATH=/root/miniconda2/bin/:/usr/bin/:$PATH

#/share/data2/guorc/Software/conda/Var/bin/samtools faidx $db
#/share/software/anaconda3/bin/bowtie2-build $db $db

if [ $Type == "Map" -o $Type == "All" ]; then

echo -e '\n\nbowtie2 & samtools running######################################################################################'
date

mkdir -p $outdir
#/share/software/anaconda3/bin/bwa mem -t $threads $db $fq1 $fq2 > $outdir/$Samp.sam
bowtie2 -u $readPair -p $threads --no-mixed -x $db -1 $fq1 -2 $fq2 --very-sensitive --n-ceil 0,0.01 -S $outdir/$Samp.sam
samtools sort -@ $threads -o $outdir/$Samp.sort $outdir/$Samp.sam

rm $outdir/$Samp.sam 
echo -e '\n\nmap finish!'

echo -e '\n\nbcftools running################################################################################################'
bcftools mpileup --threads $threads --annotate DP4 -q 20 -Q 30 -O b -f $db $outdir/$Samp.sort -o $outdir/$Samp.bcf.gz
rm $outdir/$Samp.sort

echo -e '\n\nmap & mpileup finish!'
date
fi



if [ $Type == "SNP" -o $Type == "All" ]; then

echo -e '\n\npick out MAGs (POS with >= '$SNP_depth' depth)##########################################'
date

if [ ! -f "$outdir/$Samp.vcf.filter" ];then
bcftools filter --threads $threads -g 10 -e 'TYPE=="INDEL" || (DP4[*:0]+DP4[*:1]+DP4[*:2]+DP4[*:3] < '$SNP_depth')' $outdir/$Samp.bcf.gz |\
    bcftools annotate -x INFO | bcftools view -H | awk -F ''$lab'' '{print $1"\t"$0}' > $outdir/$Samp.vcf.filter
fi

echo -e '\n\npN/pS running####################################################################################################'
date
set -e
mkdir -p $outdir/$Samp
Command="/share/data2/guorc/script/WGS/Variation/SNP_pN_pS_calculate.v7.py -f $gffPath -b $outdir/$Samp.vcf.filter -o $outdir/$Samp -t $threads -d $SNP_depth"
echo $Command
echo $Command |bash
rm $outdir/$Samp.vcf.filter


echo -e '\n\nSNP finish!'
date
fi
