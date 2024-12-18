#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-04-01, 17:15:06
# modified date: 2024-04-01, 17:15:06

'''
Jinxin Meng, 20240401
Gene Set Enrichment Analysis
Required files:
1. gene and log2foldchange of gene "\t" delimited
K00001  1.345 
K00002  2.111
K00003  -1.333
...
2. pathway_library.tsv
ko00001 K00001
ko00002 K00002
ko00003 K00004
...
'''

import sys
import pandas as pd
import blitzgsea as blitz

if len(sys.argv) != 4:
    print("Usage: calcu_gsea.py [in_f signature.tsv] [in_f pathway_library.tsv] [out_f]")
    sys.exit()

def main(signature, pathway_library, out_f):
    library = {}
    with open(pathway_library, "r") as f:
        for i in f:
            l = i.strip().split("\t")
            if l[0] not in library:
                library[l[0]] = []
            library[l[0]].append(l[1])

    signature = pd.read_csv(signature, sep = "\t")
    sorted_signature = signature.sort_values(by="1").reset_index(drop=True)
    result = blitz.gsea(sorted_signature, library)
    result.to_csv(out_f, sep="\t")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
