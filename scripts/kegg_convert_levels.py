#!/share/data1/software/miniconda3/envs/checkm2/bin/python
# -*- encoding: utf-8 -*-
##########################################################
# Creater       :  夜下凝月
# Created  date :  2020-06-26, 18:44:41
# Modiffed date :  2023-10-17, 20:39:42
##########################################################

'''
输入kegg矩阵，输出不同等级矩阵
如果一个KO被分配到多个分组，会先求平均后，再相加
如果启用--lineage T 就会输出该等级对应的上级

如：level C :      A|B|C
    levle B :      A|B

'''
import argparse
import re

PATH = "/share/data1/database/KEGG"

KEG_DICT_LEVEL_A = {}
KEG_DICT_LEVEL_B = {}
KEG_DICT_LEVEL_C = {}
KEG_DICT_LEVEL_D = {}
KEG_DICT_LEVEL_D2C = {'unknown': {'unknown'}}
KEG_DICT_LEVEL_C2B = {'unknown': {'unknown'}}
KEG_DICT_LEVEL_B2A = {'unknown': {'unknown'}}
RESULT_DICT_LEVEL_A = {}
RESULT_DICT_LEVEL_B = {}
RESULT_DICT_LEVEL_C = {}
RESULT_DICT_LEVEL_D = {}

def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", required=True, metavar='OTU', help='kegg OTU')
    parser.add_argument("-mode",  choices=['pathway', 'module'], default='pathway',
            help='convert OTU to pawhtay/module. (\033[31mdefault: pathway\033[0m)')
    parser.add_argument("--keg", required=False, default=PATH,
            help="the keg file path that inclued ko00001.keg and ko00002.keg\n\033[31mdefault: {0}\033[0m".format(PATH))
    parser.add_argument("--lineage", default="F", choices=['T','F'], help='like this: A|B|C \033[31m[T/F (default: F)]\033[0m')
    parser.add_argument("--rounding", action="store_true", help='round the value') # Jinxin Meng
    parser.add_argument("-o", required=True, metavar="output", help="output file. ")
    args = parser.parse_args()
    return args

def add_element_to_dict(dict_: dict, ko_num:str, value: str):
    """将元素添加至字典中"""
    if dict_.get(ko_num):
        dict_[ko_num].update({value})
        return True
    dict_[ko_num] = {value}

def paser_02keg(file_):
    f = open(file_, 'r')
    for line in f:
        if re.match(r"A<b>(.*)</b>", line.strip()):
            A = re.match(r"A<b>(.*)</b>", line.strip()).group(1)
        elif re.match(r"B\s+<b>(.*)</b>", line.strip()):
            B = re.match(r"B\s+<b>(.*)</b>", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_B2A, B, A)
        elif re.match(r"C\s+(.*)", line.strip()):
            C = re.match(r"C\s+(.*)", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_C2B, C, B)
        elif re.match(r"D\s+M\d{5}.*", line.strip()):
            try:
                D = re.match(r"D\s+(M\d{5})\s+(.*)\s+\[.*\]", line.strip()).group(1)
            except:
                D = re.match(r"D\s+(M\d{5})", line.strip()).group(1)
        elif re.match(r"E\s+K\d+", line.strip()):
            K = re.match(r"E\s+(K\d{5}).*", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_A, K, A)
            add_element_to_dict(KEG_DICT_LEVEL_B, K, B)
            add_element_to_dict(KEG_DICT_LEVEL_C, K, C)
            add_element_to_dict(KEG_DICT_LEVEL_D, K, D)
            add_element_to_dict(KEG_DICT_LEVEL_D2C, D, C)
    f.close()

def paser_01keg(file_):
    f = open(file_, 'r')
    for line in f:
        if re.match(r"(A\d{5})\s+(.*)", line.strip()):
            A = re.match(r"(A\d{5})", line.strip()).group(1)
        elif re.match(r"B\s+\d{5}.*", line.strip()):
            B = "B" + re.match(r"B\s+(\d{5})\s+(.*)", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_B2A, B, A)
        elif re.match(r"C\s+(\d{5})", line.strip()):
            C = "C" + re.match(r"C\s+(\d{5})\s+(.*)", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_C2B, C, B)
        elif re.match(r"D\s+K\d{5}", line.strip()):
            K = re.match(r"D\s+(K\d{5})", line.strip()).group(1)
            add_element_to_dict(KEG_DICT_LEVEL_A, K, A)
            add_element_to_dict(KEG_DICT_LEVEL_B, K, B)
            add_element_to_dict(KEG_DICT_LEVEL_C, K, C)
    f.close()

def convert(ko_name, ko_num, query_dict, result_dict, lines):
    """将ko号，转为对应的层级，并实现合并"""
    if query_dict.get(ko_name):
        names = query_dict[ko_name]
    else:
        names = {'unknown'}
    nlevels = len(names)
    for name in names:
        if result_dict.get(name):
            if result_dict[name].__len__() != ko_num.__len__():
                print(f"The ko: {ko_name}'s at \033[0,31,40m{lines}\033[m line is not equal other line's length")
                exit(0)
            result_dict[name] = list(map(lambda x, y: x + y / nlevels, result_dict[name], ko_num)) # 加和
            continue
        result_dict[name] = [ x / nlevels for x in ko_num ]

def save_result(file_, result_dict, level, tittle, lineage, rounding):
    with open(f"{file_}.level{level}", 'w', encoding='utf-8') as f:
        f.write(tittle) # 因为之前没有去除tittle的\n，所以此处不需再写\n
        for k, v in result_dict.items():
            if rounding:
                v = [str(round(x)) for x in v]
            else:
                v = [str(round(x, 4)) for x in v]
            v = "\t".join(v)
            if lineage == 'T':
                if level == "D":
                    C = list(KEG_DICT_LEVEL_D2C[k])[0]
                    B = list(KEG_DICT_LEVEL_C2B[C])[0]
                    A = list(KEG_DICT_LEVEL_B2A[B])[0]
                    k = "|".join([A, B, C, k])
                elif level == "C":
                    B = list(KEG_DICT_LEVEL_C2B[k])[0]
                    A = list(KEG_DICT_LEVEL_B2A[B])[0]
                    k = "|".join([A, B, k])
                elif level == "B":
                    A = list(KEG_DICT_LEVEL_B2A[k])[0]
                    k = "|".join([A, k])
            f.write(f"{k}\t{v}\n")

def main(OTU_file, output_file, mode, lineage, rounding): # Jinxin Meng
    f = open(OTU_file, 'r')
    lines = 1
    print(f"process and calculate ...")
    tittle = f.readline() # 跳过 tittle
    for line in f:
        line_split = re.split(r"\t", line.strip())
        ko_name = line_split[0]
        #ko_num = [int(x) for x in line_split[1:]]
        ko_num = [float(x) for x in line_split[1:]]
        convert(ko_name, ko_num, KEG_DICT_LEVEL_C, RESULT_DICT_LEVEL_C, lines)
        convert(ko_name, ko_num, KEG_DICT_LEVEL_B, RESULT_DICT_LEVEL_B, lines)
        convert(ko_name, ko_num, KEG_DICT_LEVEL_A, RESULT_DICT_LEVEL_A, lines)
        if mode == "module":
            convert(ko_name, ko_num, KEG_DICT_LEVEL_D, RESULT_DICT_LEVEL_D, lines)
        lines += 1
    f.close()
    print(f"save result ...")
    save_result(output_file, RESULT_DICT_LEVEL_A, 'A', tittle, lineage, rounding) # Jinxin Meng
    save_result(output_file, RESULT_DICT_LEVEL_B, 'B', tittle, lineage, rounding) # Jinxin Meng
    save_result(output_file, RESULT_DICT_LEVEL_C, 'C', tittle, lineage, rounding) # Jinxin Meng
    if mode == "module":
        save_result(output_file, RESULT_DICT_LEVEL_D, 'D', tittle, lineage, rounding) # Jinxin Meng


if __name__ == "__main__":
    args = get_args()
    input_file = args.i
    output_file = args.o
    path = args.keg
    mode = args.mode
    lineage = args.lineage

    if mode == 'pathway':
        paser_01keg(f"{path}/ko00001.keg")
    elif mode == 'module':
        paser_02keg(f"{path}/ko00002.keg")
    else:
        print("not the keg {mode}")
        exit(127)
    main(input_file, output_file, mode, lineage, args.rounding) # Jinxin Meng
