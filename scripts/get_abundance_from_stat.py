#!/share/data1/zhangy/software/miniconda/envs/bio37/bin/python
# -*- encoding: utf-8 -*-
import pandas as pd
import argparse

def calc_TPM(df, output_file) :
    df['RPK'] = df['mapping'].div(df['length']/1000)
    temp_sum = df['mapping'].div(df['length']/1000).sum()/(10**6)
    df['TPM'] = df['RPK'].div(temp_sum)
    df.to_csv(output_file, header=True, index=False, sep="\t")
    print("Finist:", output_file)

def calc_RPKM(df, output_file):
    total_reads = df['mapping'].sum()+df['unmapping'].sum()
    df['RPKM'] = df['mapping']/((df['length']-150)*total_reads)*(10**9)
    df.to_csv(output_file, header=True, index=False, sep="\t")
    print("Finish:",output_file)

def main(input_file=None, output_file=None):
    '''计算RPKM'''
    df = pd.read_csv(input_file, header=None, sep="\t")
    dict_name = {0:'name', 1:'length', 2:'mapping', 3:'unmapping'}
    df = df.rename(columns=dict_name)

    calc_TPM(df, output_file)
    #  calc_RPKM(df, output_file)

if __name__ == "__main__":
    '''设置帮助信息'''
    parser = argparse.ArgumentParser(description="calculate abundance acording to stat file")
    parser.add_argument("i",help="input_stat_file")
    parser.add_argument("-o",help="output_file [default: source file]")
    args = parser.parse_args()

    input_file = args.i
    output_file = args.o

    if args.o == None:
        output_file = args.i

    main(input_file=input_file, output_file=output_file)
