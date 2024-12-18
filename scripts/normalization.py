#!/share/data1/software/miniconda3/bin/python3
import pandas as pd
import numpy as np
import sys

"""
归一化profile， 抽平样本reads数目,抽平的reads数由脚本自己取最小值
表中数据将被转为整数
"""
def check_args():
    if sys.argv.__len__() > 5 or sys.argv.__len__() < 4 or '-h' in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
        print(f"{sys.argv[0]}  <profile>  <output> <seed> [minimum: default is auto]")
        print("\033[31m<output> wile be overwritten\033[0m")
        exit(192)

def process():
    stat=1

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    seed = int(sys.argv[3])
    args_len = len(sys.argv)

    np.random.seed(seed)

    #df = pd.read_csv(in_file, sep="\s+", index_col=0, header=0).astype(int)
    df = pd.read_csv(in_file, sep="\t", index_col=0, header=0).astype(int)

    try:
        Min_num = int(sys.argv[4])
    except:
        print("You do not set the minimum value, it will be discovered by auto.")
        Min_num = df.sum().min()

    col_names = df.columns
    INDEX = df.index

    result_list = []

    for i in col_names[0:]:
        print(f"process sample {stat}: {i}")
        stat+=1
        if args_len == 4 or df[i].sum() >= Min_num:
            temp = np.array(df[i])
            temp1 = np.random.choice(np.repeat(INDEX, temp),Min_num)
            temp1 = pd.DataFrame(temp1)[0].value_counts()
            temp1 = pd.DataFrame(temp1)
            temp1.columns = [i]
        else:
            temp1 = pd.DataFrame(df[i])
        result_list.append(temp1)

    result = pd.concat(result_list, join='outer', axis=1).fillna(0).astype(int)
    result.to_csv(out_file, sep="\t")
    print(f"min reads is {Min_num}")

if __name__ == "__main__":
        check_args()
        process()
