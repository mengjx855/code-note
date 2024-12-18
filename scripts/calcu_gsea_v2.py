#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-04-01, 17:15:06
# modified date: 2024-04-01, 18:28:59

'''
Jinxin Meng, 20240401
Gene Set Enrichment Analysis
Required files:
1. profile
name    s1  s2  ...
K00001  392 230
K00002  122 340
...
2. sample_group
sample  group
s1      x
s2      y
3. group_pair
x   y
x   z
4. pathway_library.tsv
ko00001 K00001
ko00002 K00002
ko00003 K00004
...
'''

import sys
import pandas as pd
import numpy as np
import blitzgsea as blitz

if len(sys.argv) != 6:
    print("Usage: calcu_gsea_v2.py [in_f profile] [in_f sample_group] [in_f group_pair] [in_f pathway_library.tsv] [out_prefix]")
    sys.exit()

def parse_library(pathway_library):
    library = {}
    with open(pathway_library, "r") as f:
        for i in f:
            l = i.strip().split("\t")
            if l[0] not in library:
                library[l[0]] = []
            library[l[0]].append(l[1])
    return library

def parse_gp(group_pair):
    gp = []
    with open(group_pair, "r") as f:
        for i in f:
            l = i.strip().split("\t")
            gp.append(l)
    return gp

def signature_rank(profile, group, gp):
    name = list(profile.index)
    sample_x = list(group.loc[group["group"] == gp[0], "sample"])
    sample_y = list(group.loc[group["group"] == gp[1], "sample"])
    mean_x = profile.loc[:, sample_x].mean(axis = 1) + 0.0001
    mean_y = profile.loc[:, sample_y].mean(axis = 1) + 0.0001
    res = pd.DataFrame({"0": name, "1": list(np.log2(mean_x / mean_y))}).sort_values(by="1", ascending=False).reset_index(drop=True)
    return res

def main(profile, group, group_pair, pathway_library, out_prefix):
    library = parse_library(pathway_library)
    profile = pd.read_csv(profile, sep="\t", header = 0, index_col = 0)
    group = pd.read_csv(group, sep="\t", header = 0)
    gp = parse_gp(group_pair)
    
    for i in gp:
        signature = signature_rank(profile, group, gp=i)
        result = blitz.gsea(signature, library)
        result.to_csv(out_prefix + "_" + i[0] + "_vs_" +  i[1] + ".tsv", sep="\t")
        fig_table = blitz.plot.top_table(signature, library, result, n=150)
        fig_table.savefig(out_prefix + "_" + i[0] + "_vs_" +  i[1] + ".figtable.pdf")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
