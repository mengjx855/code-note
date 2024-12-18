#!/share/data1/zhangy/software/miniconda3/envs/bio37/bin/python
# -*- encoding: utf-8 -*-

###########################################################
# Author:				ZhangYue
# Description:				
# Time:				2019年07月30日	 Tuesday
###########################################################
import re
import os
import sys
import argparse
def main(input_path=None,output_path=None,database=None):
    samtools = "/share/data1/zhangy/software/miniconda3/envs/bio37/bin/samtools"
    bowtie2 = "/share/data1/zhangy/software/miniconda3/envs/bio37/bin/bowtie2"
    sh_content =""
    info_dict = {}
    info_dict['output_path'] = output_path
    info_dict['database'] = database
    info_dict['samtools'] = samtools
    info_dict['bowtie2'] = bowtie2
    bam_sort_stat = ""
    j = 1
    f = open(input_path)
    for i in f:
        i = i.strip()
        temp = re.split("\s+", i)
        info_dict['sample_name'] = temp[0]
        info_dict['fq1'] = temp[-2]
        info_dict['fq2'] = temp[-1]
        sh_content += "date\n"
        sh_content += "mkdir -p {output_path}/{sample_name}\n".format(**info_dict)
        sh_content += "cd {output_path}/{sample_name}\n".format(**info_dict)
        sh_content += "{bowtie2} --end-to-end --very-fast --threads 20 -x {database} -1 {fq1} -2 {fq2} --al-conc {sample_name}.al.con.gz -S {sample_name}.sam 2> {sample_name}.log\n".format(**info_dict)
        sh_content += "{samtools} view -@ 5 -bS {sample_name}.sam | {samtools} sort - -m 10G -o {sample_name}.bam.sort\n".format(**info_dict)# sam -> bam
        sh_content += "ls\n"
        sh_content += "{samtools} index {sample_name}.bam.sort\n".format(**info_dict)# bam.sort -> bam.sort.stat
        sh_content += "ls\nls\nls\n"
        sh_content += "{samtools} idxstats {sample_name}.bam.sort > {sample_name}.bam.sort.stat\n".format(**info_dict)# bam.sort -> bam.sort.stat
        sh_content += 'if [ -s "{sample_name}.bam.sort.stat" ];then echo "rm {sample_name}.sam";rm {sample_name}.sam ;fi\n\n\n'.format(**info_dict)
        bam_sort_stat += "{sample_name}\t{output_path}/{sample_name}/{sample_name}.bam.sort.stat\n".format(**info_dict)
# -------------------------------------^90-----------------------------------------------------------
# -------------------------------------80------------------------------------------------------------

        #  sh_content += "/share/data7/zhangy2/Software/.conda/envs/bio37/bin/bowtie2 --end-to-end --sensitive --score-min L,-1.2,-1.2 --threads 20 -x {database} -1 {fq1} -2 {fq2} --al-conc {sample_name}.al.con.gz -S {sample_name}.sam 2> {sample_name}.log\n".format(**info_dict)
        #  sh_content += "/share/data7/zhangy2/Software/.conda/envs/bio37/bin/samtools view -@ 5 -bS {sample_name}.sam |  /share/data7/zhangy2/Software/.conda/envs/bio37/bin/samtools sort - -m 10G -o {sample_name}.bam.sort\n".format(**info_dict)# sam -> bam
        #  sh_content += "ls\n"
        #  sh_content += "/share/data7/zhangy2/Software/.conda/envs/bio37/bin/samtools index {sample_name}.bam.sort\n".format(**info_dict)# bam.sort -> bam.sort.stat
        #  sh_content += "ls\nls\nls\n"
        #  sh_content += "/share/data7/zhangy2/Software/.conda/envs/bio37/bin/samtools idxstats {sample_name}.bam.sort > {sample_name}.bam.sort.stat\n".format(**info_dict)# bam.sort -> bam.sort.stat
        #  sh_content += "date\n"
        #  sh_content += 'if [ -s "{sample_name}.bam.sort.stat" ];then echo "rm {sample_name}.sam";rm {sample_name}.sam ;fi\n\n\n'.format(**info_dict)
        print(j,":",i)
        j += 1
    with open("{output_path}/mapping_stat.list".format(**info_dict), 'w', encoding='utf-8') as f:
        f.write(bam_sort_stat)
    with open("{output_path}/mapping_and_stat.sh".format(**info_dict),'w',encoding="utf-8") as f:
        f.write(sh_content)
        print("Finish")
        print("{output_path}".format(**info_dict))
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-db', metavar='database', type=str, help='the database index')
    parser.add_argument('-i', metavar='if', type=str, help='the file include: sample \\t fq1 \\t fq2')
    parser.add_argument('-o', metavar='of', type=str, help='output_path')
    if sys.argv.__len__() != 7:
        parser.print_help()
        exit(0)
    #  args = parser.parse_args("-h".split())
    args = parser.parse_args()

    input_path = args.i
    output_path = os.path.abspath(args.o)
    database = os.path.abspath(args.db)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    main(input_path=input_path,output_path=output_path,database=database)




