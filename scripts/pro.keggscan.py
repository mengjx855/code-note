#!/share/data1/zhangy/software/miniconda3/envs/bio37/bin/python
# -*- encoding: utf-8 -*-
###########################################################
# Author:				ZhangYue
# Description:				
# Time:				2020年05月14日	 Thursday
###########################################################
import re
import sys


def main(input_file_list,  min_score):
    global RESULT
    f = open(input_file_list, 'r')
    # 跳过前两行
    f.readline()
    f.readline()
    last_gene = ""
    last_line = ""
    last_score = ""

    for line in f:
        line_split = re.split(r"\s+", line.strip())
        if line[0] == "*":
            gene, _, _, score, _ = line_split[1:6]
            line = re.match(r"\*\s+(.*)", line.strip()).group(1)
        else:
            gene, _, _, score, _ = line_split[0:5]

        score = float(score)

        if score < min_score:
            continue

        if last_gene != gene:
            if last_gene != "":
                print(last_line)
            last_gene = gene
            last_line = line.strip()
            last_score = score
        elif last_gene == gene:
            if last_score < score:
                last_line = line.strip()
                last_score = score
    print(last_line)
    f.close()



if __name__ == "__main__":
    RESULT = {}
    if sys.argv.__len__() != 3:
        print(f"{__file__} <keggscan_result>  <score>")
        exit(0)
    main(sys.argv[1], float(sys.argv[2]))
