# kegg
# ko00001.keg | pathway 的从属关系 - Last updated: April 3, 2023
# ko00002.keg | module 的从属关系 - Last updated: March 29, 2023
# br08610.keg | kegg物种信息和ncbi的对应关系 - Last updated: April 3, 2023
# org_list.html | organisms kegg物种信息 - Last updated: April 3, 2023

wget -O "org_list.html" "https://www.kegg.jp/kegg/catalog/org_list.html"
wget -O "ko00001.keg" "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir="
wget -O "ko00002.keg" "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
wget -O "br08610.keg" "https://www.kegg.jp/kegg-bin/download_htext?htext=br08610.keg&format=htext&filedir="
# 这个URL非常有用 https://www.kegg.jp/kegg-bin/download_htext?htext=br08610.keg&format=htext&filedir=
# 这个连接就是存储各类层级文件的一个链接，*keg就是文件格式，找到这个文件就可以进行下载
# 例如，我现在检索mmu05416通路，返回的网页是中点击pathway menu，地址栏中出现https://www.genome.jp/brite/query=05416&htext=br08901.keg&option=-a&node_proc=br08901_org&proc_enabled=mmu&panel=collapse
# 其中，br08901.keg就是存储小鼠通路层级关系的文件，根据上述链接进行下载和处理
# 当然，这些文件的编号可以在https://www.genome.jp/kegg/brite.html这里找到
le br08901.keg | perl -lne 'if(/^A<b>(.*?)</){$a=$1}elsif(/^B\s+(.*)/){$b=$1;}elsif(/^C\s+(\d+)\s+(.*)$/){print"$a\t$b\tmap$1\t$2"}' | csvtk.exe add-header -t -n "lvA,lvB,lvC,lvCdes" > br08901.tsv2

# 处理br08610.keg
scripts/get_taxa_id.pl br08610.keg kegg_taxa_id
# 处理org_list.html
scripts/get_taxa_info.pl org_list.html kegg_taxa_info
# 处理生成下载命令
scripts/generate_dwn_script.pl kegg_taxa_info r2.dwn_pep.sh r3.dwn_keg.sh

# 下载数据
# 有的数据下载不了，但是尽量多下载
mkdir kegg_KO ncbi_pep
parallel -j 24 --joblog ncbi_pep_parallel_joblog.log < r2.dwn_pep.sh > ncbi_pep.log 2>&1 &
parallel -j 24 --joblog kegg_KO_parallel_joblog.log < r3.dwn_keg.sh > kegg_KO.log 2>&1 &

# 检查没下载到的数据
find ncbi_pep/ -name "*" -type f -size 0c
find kegg_KO/ -name "*" -type f -size 0c

# 匹配蛋白序列和KO 
mkdir pep_KO
cut -f1 kegg_taxa_info | awk '{print "gzip -d ncbi_pep/"$0".pep.fasta.gz ; python scripts/process_proteins.py --org "$0" --pep ncbi_pep/"$0".pep.fasta --keg kegg_KO/"$0"00001.keg --out pep_KO/"$0"_KO.fasta ; gzip ncbi_pep/"$0".pep.fasta"}' > r4.match.sh
parallel -j 32 --joblog pep_KO_parallel_joblog.log < r4.match.sh > pep_KO.log 2>&1 &
find pep_KO/ -name "*" -type f -size 0c | xargs rm

# 下载module信息 | r5.download_module.sh
./scripts/parse_module.pl ko00002.keg > ko00002.txt
mkdir module_map module_KO
awk -F '\t' '{print $4}' ko00002.txt | perl -ne 'chomp;print "wget -c --timeout 60 --tries 10 -O \"module_map/$_.html\" \"https://www.genome.jp/dbget-bin/www_bget?md:$_\" \|\| rm \"module_map/$_.html\"\n"' > r5.dwn_module_map.sh
awk -F '\t' '{print $4}' ko00002.txt | perl -ne 'chomp;print "wget -c --timeout 60 --tries 10 -O \"module_KO/$_.html\" \"https://www.kegg.jp/kegg-bin/module_ko_list?map=$_&org=ko\" \|\| rm \"module_KO/$_.html\"\n"' > r6.dwn_module_KO.sh
parallel -j 24 --joblog module_map_parallel_joblog.log < r5.dwn_module_map.sh > module_map.log 2>&1 &
parallel -j 24 --joblog module_KO_parallel_joblog.log < r6.dwn_module_KO.sh > module_KO.log 2>&1 &


# 生成maplist
ls --color=auto module_map/*.html | parallel --plus -j 32 sh ./scripts/modules.map.parse.sh {} {/.} > module_map.list
#替换错误字符串
sed -i 's/<wbr>//g' module_map.list
sed -i 's/ --//g' module_map.list
sed -i 's/-- //g' module_map.list
# 提取反应式
 erl scripts/kegg_module_breaker.pl module_map.list module_map.reaction

# 提取ko对应信息
le ko00001.keg | perl -ne 'chomp; next if !/^D\s+(K\d+)/;$_=~/^D\s+(K\d+)\s+(.*)/;print "$1\t$2\n"' | sort -k 1,1 -u > KO_description

# 合并数据
cat pep_KO/* > KEGG20230401.fas

# diamond 建库
diamond makedb --threads 40 --in KEGG20230401.fas -d KEGG20230401.dmnd

./scripts/parse_KO.pl ko00001.keg > KO_level_A_B_C_D_Description
./scripts/parse_module.pl ko00002.keg > module_level_A_B_C_D_Description

ls *html | cut -d. -f1 | parallel -j 10 -q perl -e 'open I,$ARGV[1];while(<I>){chomp;if(/>(K\d+)<\/a/){print "$ARGV[0]\t$1\n"}}' {} {}.html > ../module_KO.list &