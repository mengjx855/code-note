#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-03-12, 11:55:10
# modified date: 2024-03-12, 18:13:31

'''
Jinxin Meng, 20240312
Generate the adjacent matrix to determinethe presence or absence of feature in each sample.
All genes in each gene cluster were regarded have same feature annotation.
Like KO, CAZyme, COGs ... annotation were used in this script.
Example out:
name    s1  s2  s3 ...
K00001  1   4   2       
K00002  2   1   1
K00003  1   2   3
...
Required files:
1. KO annotation with field gene name and feature name with "\t" or "\s" delimited
2. gene cluster relationship provided by mmseqs
'''

import sys

if len(sys.argv) != 4:
    print("Usage: get_adjacency_for_samples.py [in_f <gene ID|feature ID>] [gene_cluster.tsv <rep. gene ID|gene ID>] [out_f]")
    print("  Reqired files:")
    print("  1. in_f with field gene ID and feature ID, and gene ID with reg. exp. \"(.*?_\S+)\", such as s1_K141_100_1")
    print("  2. gene_cluster.tsv with field representative gene ID and all gene ID ")
    sys.exit()

def main(in_f, cls_f, out_f):
    feature = {}
    with open(in_f, "r") as f:
        for i in f:
            l = i.strip().split()
            feature[l[0]] = l[1]
    
    res = {}
    sample = []
    with open(cls_f, "r") as f:
        for i in f:
            l = i.strip().split()
            if l[0] not in feature:
                continue
            x = l[1].split("_")[0]
            if feature[l[0]] not in res:
                res[feature[l[0]]] = {}
            if x not in res[feature[l[0]]]:
                res[feature[l[0]]][x] = 0
            if x not in sample:
                sample.append(x)
            res[feature[l[0]]][x] += 1
    
    sample = sorted(sample)
    out_f = open(out_f, "w")
    out_f.write("name\t" + "\t".join(sample) + "\n")
    for k, v in res.items():
        l = [k]
        for i in sample:
            if i in v:
                l.append(str(v[i]))
            else:
                l.append("0")
        out_f.write("\t".join(l) + "\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])