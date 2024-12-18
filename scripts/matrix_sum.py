#!/share/data1/zhangy/software/miniconda/envs/bio37/bin/python
# -*- encoding: utf-8 -*-
###########################################################
# Author:				ZhangYue
# Description:				See below
# Time:				2020年07月06日	 Monday
# ModiffTime:       2020年09月01日
###########################################################
'''
Version : V2.0

group: Must non header
    name_1  group_1     Group_A
    name_1  group_4     Group_B
    name_2  group_2     Group_C
    name_3  group_3     Group_A
    ...     ...         ...

matrix: Must have header
    tittle sample_1 sample_2    ...
    name_1  1       1           0
    name_2  0       1           2
    name_3  3       4           9
    ...     ...     ...         ...
'''
import  argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="Matrix", required=True, help="Input matrix")
    parser.add_argument("-g", metavar="group file", required=True, help="Input group")
    parser.add_argument("-V", metavar="Var col", default = 1, type=int, help="The ID's columns in group [default: 1]")
    parser.add_argument("-G", metavar="Group col", default = 2, type=int, help="The group's columns in group file [default: 2]")
    parser.add_argument("-o", metavar="output", required=True, help="Output matrix")
    parser.add_argument("-s", required=False, type=int,  default=1, choices=[0,1,2,3], help="Split Str fo matrix and group.[default: 1] \033[31mThey must have similary split str\033[0m\n0 -> \\t\n1 -> \\s+\n2 -> |\n3 -> ,")
    args = parser.parse_args()
    return args

def parse_group(file_, GRP, sstr, V, G):
    '''
    只是解析分组文件，{key: {taxo_1, taxo_2,...}, ...}
    '''
    f = open(file_, 'r')
    for line in f:
        line_split = re.split(sstr, line.strip())
        if GRP.get(line_split[V]):
            GRP[line_split[V]].update({line_split[G]})
            continue
        GRP[line_split[V]] = {line_split[G]}
    f.close()

def main(file_, sstr, output_file):
    result = {}
    f = open(file_, 'r')
    tittle = f.readline()
    for line in f:
        line_split = re.split(sstr, line.strip())
        name = line_split[0]
        try:
            groups = GRP[name] # 将名字转为分组 这边得到一个set，后面需要循环才行
        except:
            print(f"\033[31m{name}\033[0m can't find in group file, it will be changed to unknown.")
            groups = {'unknown'}
        num = [float(x) for x in line_split[1:]] # 如果需要整数，此处改为int即可
        for group in groups:
            if result.get(group):
                result[group] = list(map(lambda x,y: x+y, result[group], num))
            else:
                result[group] = num
    f.close()

    with open(output_file, 'w', encoding="utf-8") as f:
        f.write(tittle)
        if sstr != "\|":
            sstr = "\t"
        for k,v in result.items():
            v = [str(x) for x in v]
            f.write(f"{k}{sstr}{sstr.join(v)}\n")


if __name__ == "__main__":
    args = get_args()
    GRP = {}
    sstr = {0:"\t", 1:"\s+", 2:"\|", 3:","}[args.s]
    if args.V < 0 or args.G < 0 or args.V == args.G:
        print("Plese check Your parameter -V -G, they must >= 1 and not equal!!!")
        exit(127)
    V = args.V - 1 
    G = args.G - 1 
    parse_group(args.g, GRP, sstr, V, G)
    main(args.i, sstr, args.o)
