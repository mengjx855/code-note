#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-12-21, 22:04:06
# modified date: 2023-12-22, 21:19:27

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-bg', type=str, metavar='background', help='gff-like table with field contigs|gene|start|end|direction.')
    parser.add_argument('-g1', type=str, metavar='ARGs_gene_list', help='a table containing target_gene list in first column.')
    parser.add_argument('-g2', type=str, metavar='MGEs_gene_list', help='a table containing target_gene list in first column.')
    parser.add_argument('-l', type=int, metavar='interval_len', default=5000, help='interval lenth for up- and down-stream of target gene')
    parser.add_argument('-o', type=str, metavar='out_prefix', help='output graph file')
    args = parser.parse_args()
    return args

def get_arg_gene(in_f):
    gene_dict = {}
    with open (in_f, 'r') as file:
        for l in file:
            gene = l.split()[0]
            contig = re.findall(r'(\S+)_\d+$', gene)[0]
            if gene not in gene_dict:
                gene_dict[gene] = contig
            else:
                continue
    return gene_dict

def get_mge_gene(in_f):
    gene_list = []
    with open (in_f, 'r') as file:
        for l in file:
            gene = l.split()[0]
            if gene not in gene_list:
                gene_list.append(gene)
            else:
                continue
    return gene_list


def get_bg_info(in_f, ref_contig_list):
    bg_gene_dict = {}
    bg_interval_dict = {}
    with open (in_f, 'r') as file:
        for l in file:
            contig, gene, start, end, direction = l.split()
            if contig in ref_contig_list:
                if contig not in bg_gene_dict:
                    bg_gene_dict[contig] = []
                    bg_interval_dict[contig] = []
                bg_gene_dict[contig].append(gene)
                bg_interval_dict[contig].append([int(start), int(end), str(direction)])
            else:
                continue
    return bg_gene_dict, bg_interval_dict

def out_g(out_dict, out_g, mge_gene_list):
    out_g = open(out_g, 'w')
    for i in out_dict.keys():
        seqs = sorted(out_dict[i])
        out_g.write(i + '\t')
        for j in seqs:
            x = j
            if x == i:
                x += ' (ARGs)'
            elif x in mge_gene_list:
                x += ' (MGEs)'
            if out_dict[i][j][2] == '+':
                out_g.write('|== ' + x + ' ==>|')
            else:
                out_g.write('|<== ' + x + ' ==|')
        out_g.write('\n')
    out_g.close()
    

def out_f(out_dict, out_f, mge_gene_list):
    out_f = open(out_f, 'w')
    for i in out_dict.keys():
        seqs = sorted(out_dict[i])
        for j in seqs:
            x = j
            if x == i:
                x += '\tARGs'
            elif x in mge_gene_list:
                x += '\tMGEs'
            else:
                x += '\tOther'
            out_f.write('\t'.join([i, x]) + '\t' + '\t'.join(str(s) for s in out_dict[i][j]) + '\n')
    out_f.close()

def main(bg, g1, g2, l, o):
    # 获得耐药基因基因和contig的信息 gene_dict[gene] = contig
    gene_dict = get_arg_gene(g1)
    contig_list = list(set(gene_dict.values()))
    # MGE基因列表
    mge_gene_list = get_mge_gene(g2)
    # 背景基因信息 bg_gene_dict[contig] = [gene_name] bg_interval_dict[contig] = [gene_loc]
    bg_gene_dict, bg_interval_dict = get_bg_info(bg, contig_list)
    
    # 滑动基因进行ARGs上下游某个区间判断基因
    out_dict = {}
    for k, v in gene_dict.items(): 
        gene_index = bg_gene_dict[v].index(k)
        gene_interval = bg_interval_dict[v][gene_index]
        n_gene = len(bg_gene_dict[v])
        
        '''
        1. 如果目的基因单独存在某个contig上，不进行区间判断
        2. 如果目的基因位于某个contig的最左端，只判断另一侧 
        3. 如果目的基因位于某个contig的最右端，只判断另一侧
        4. 如果目的基因位于contig的中间，两侧均判断
        '''
        out_dict[k] = {}
        out_dict[k][k] = gene_interval
        i = 1
        if gene_index == 0 and n_gene == 1:
            continue
        elif gene_index == 0 and n_gene > 1:
            lim_right = gene_interval[1] + l
            while bg_interval_dict[v][gene_index + i][0] <= lim_right:
                out_dict[k][bg_gene_dict[v][gene_index + i]] = bg_interval_dict[v][gene_index + i]
                i+=1
                if gene_index + i > (n_gene - 1):
                    break
        elif gene_index == (n_gene - 1) and n_gene > 1:
            lim_left = gene_interval[0] - l
            while bg_interval_dict[v][gene_index - i][1] >= lim_left:
                out_dict[k][bg_gene_dict[v][gene_index - i]] = bg_interval_dict[v][gene_index - i]
                i+=1
                if gene_index - i < 0:
                    break
        else:
            lim_left = gene_interval[0] - l
            lim_right = gene_interval[1] + l
            while bg_interval_dict[v][gene_index - i][1] >= lim_left:
                out_dict[k][bg_gene_dict[v][gene_index - i]] = bg_interval_dict[v][gene_index - i]
                i+=1
                if gene_index - i < 0:
                    break
            i = 1
            while bg_interval_dict[v][gene_index + i][0] <= lim_right:
                out_dict[k][bg_gene_dict[v][gene_index + i]] = bg_interval_dict[v][gene_index + i]
                i+=1
                if gene_index + i > (n_gene - 1):
                    break

    # 判断MGE基因是否在这些范围内，在的话，就保留下来这个out_dict的子字典
    # 做集合运算
    # 循环修改字典时，迭代对象不能时字典
    y = set(mge_gene_list)
    for k in list(out_dict.keys()): 
        x = set(out_dict[k].keys())
        if len(x & y) == 0:
            del out_dict[k]
    
    # 输出结果
    out_f(out_dict, o + '_find.tsv', mge_gene_list)
    out_g(out_dict, o + '_graph.tsv', mge_gene_list)

if __name__ == '__main__':
    args = get_args()
    main(args.bg, args.g1, args.g2, args.l, args.o)

