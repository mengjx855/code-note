#!/share/data1/software/miniconda3/envs/jinxin/bin/python3

from scipy import stats
import numpy as np
import pandas as pd
import sys

print(f"read files: {sys.argv[1]}, {sys.argv[2]}")
df1 = pd.read_csv(f"{sys.argv[1]}", sep="\t",header=0,index_col=0)
df2 = pd.read_csv(f"{sys.argv[2]}", sep="\t", header=0, index_col=0)

# align column names
df2 = df2.reindex(columns=df1.columns)
dt1 = df1.rank(axis=1)
dt2 = df2.rank(axis=1)

nrow1 = dt1.shape[0]
nrow2 = dt2.shape[0]

out_f = f"{sys.argv[3]}"
f = open(out_f, "w")

otus1 = dt1.index
otus2 = dt2.index

f.write("name_a\tname_b\tcorr\tpvalue\n")
for i in range(0,nrow1):
    for j in range(0,nrow2):
        s = stats.pearsonr(dt1.iloc[i,:], dt2.iloc[j,:])
        temp_str = f"{otus1[i]}\t{otus2[j]}\t{s[0]}\t{s[1]}\n"
        f.write(temp_str)
f.close()
