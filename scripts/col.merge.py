#!/share/data1/zhangy/software/miniconda/envs/bio37/bin/python
# -*- encoding: utf-8 -*-
###########################################################
# Author:			YeXiaNinaGyue
# Description:		根据列存在缺失合并矩阵
# Time:				2020年11月20日	 Friday
###########################################################
import  argparse
import re
import copy

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', metavar='Matrix', type=str, help='matrix')
    parser.add_argument('-g', metavar='Group', type=str, help='Sample group')
    parser.add_argument('-o', metavar='output', type=str, help='output file')
    parser.add_argument('-p', metavar='percent', type=float, help=r'gatter than 0 is 1 for per group')
    args = parser.parse_args()
    return args

def parse_group(file_, output_dict):
    f = open(file_, 'r')
    groups = {}
    for line in f:
        ls = re.split(r"\t", line.strip())
        groups[ls[0]] = ls[1]
        output_dict[ls[1]] = "0"
    f.close()
    return groups

def judge(title, ls, temp_dict, percent, groups):
    '''
    只返回大于等于阈值的分组字典{group1: "1", group5: "1", ...}
    '''
    tt = {}
    for i, v in enumerate(ls[1:], 1):
        group = groups[title[i]]
        if v != "0" and v!= "0.0":
            if temp_dict.get(group):
                temp_dict[group][0] += 1
                temp_dict[group][1] += 1
            else:
                temp_dict[group] = [1, 1] # 第一个用于计数非0的个数， 第二个计数分组有多少样本
        else:
            if temp_dict.get(group):
                temp_dict[group][1] += 1
            else:
                temp_dict[group] = [0, 1] # 第一个用于计数非0的个数， 第二个计数分组有多少样本
    for k,v in temp_dict.items():
        if v[0]/v[1] >= percent:
            tt[k] = "1"
    return tt

def main(input_file, groups, percent, output, output_dict):
    f = open(input_file, 'r')
    title = re.split(r"\s+", f.readline().rstrip())
    title = dict(enumerate(title)) # {value: index}
    output_file = open(output, 'w')
    output_file.write("Familys\t{}\n".format("\t".join(output_dict.keys())))
    for line in f:
        temp_output_dict = copy.deepcopy(output_dict)
        temp_dict = {}
        ls = re.split("\s+", line.strip())
        temp_dict = judge(title, ls, temp_dict, percent, groups) # 根据返回值修改temp_outpu_dict
        for k,v in temp_dict.items():
            temp_output_dict[k] = v
        #output_file.write("{}\t{}\n{}\t{}\n".format(ls[0], "\t".join(temp_output_dict.values()), ls[0], "\t".join(temp_output_dict.keys())))
        output_file.write("{}\t{}\n".format(ls[0], "\t".join(temp_output_dict.values())))
    f.close()
    output_file.close()


if __name__ == "__main__":
    args = get_args()
    output_dict = {}
    groups = parse_group(args.g, output_dict)
    main(args.i, groups, args.p, args.o, output_dict)
