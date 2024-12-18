#!/usr/bin/bash
if [ $# -lt 4 ];then
    echo "$0 <out> <threads> <block_size:bytes> <url>"
    exit 127
fi

out_f=$1
threads=$2
size_split=$3 # 就是文件将分为多大的块进行逐个下载
url=$4
if [ -f $out_f ];then
    echo "${out_f} is exists."
    exit 127
fi

temp_dir=${out_f}_$$  # $$ is pid of this script 
file_size=`curl -I -s  ${url} | perl -ne 'if (/Content-Length:\s+(\d+)/){print "$1";last}'`;  #获取文件的具体大小 byte

mkdir ${temp_dir};
perl -e 'foreach ($i=0; $i<=$ARGV[0]; $i+=$ARGV[1]){$s=$i+1; $s=0 if $i == 0; $e=$i+$ARGV[1];$e=$ARGV[0] if $e > $ARGV[0];print "$s\-$e\n"}'  ${file_size} ${size_split} > ${temp_dir}/sp.list 
parallel -j ${threads} --joblog ${temp_dir}/parallel.log  curl -r {} -o ${temp_dir}/{\#}  ${url} :::: ${temp_dir}/sp.list ;

nfiles=`cat ${temp_dir}/sp.list| wc -l`;
# 如果上步出错，就重试10次，间隔5秒
count_=0
while [ $? -ne 0 ]
do
    $count+=1
    sleep 5
    parallel --retry-failed --joblog ${temp_dir}/parallel.log ;
    if [ $count -lt 10 ];then
        echo "faild"
        echo "Log file: ${temp_dir}/parallel.log"
        echo "You Can run: parallel --retry-failed --joblog ${temp_dir}/parallel.log to retry downolad"
        echo "And then run : for i in \`seq 1 $nfiles\`;do cat ${temp_dir}/$i >> ${out_f};done   to combine."
        exit 127
    fi
done
touch ${out_f};
for i in `seq 1 $nfiles`;do cat ${temp_dir}/$i >> ${out_f};done
