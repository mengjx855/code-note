#!/usr/bin/env python

import fileinput
import os, sys
import collections

prev = []
result = []

if len(sys.argv) < 2: 
    print "Usage: conn_blast.py blast_fmt6 blast_fmt6.conn"
    sys.exit()

input_file = open(sys.argv[1])
output_file = sys.argv[2]

oot = open(output_file, 'w')

## Required functions
def convert(str_a):
    try:
        str_a = float(str_a)
    except ValueError:
        str_a = str_a
    return(str_a)

def reverse_se(x):
    a = x[6]
    b = x[7]
    if  int(a) > int(b):
        a,b = b,a
        x[6] = a
        x[7] = b
    return(x)

## load table
line_list = []
for line in input_file:
    line = line.strip().split("\t")
    line = map(lambda x:convert(x), line)
    line_list.append(line)

# print line_list
## reset the start and sort lines (reverse the oeder when sstart < send)
line_list_reset = map(lambda x:reverse_se(x), line_list)
sort_lines = sorted(line_list_reset, key = lambda x:[x[0], x[1], x[6], x[7], -x[-1]])

# for i in sort_lines:
#     print i
for s in sort_lines:
    blast_id = s[2]
    bitscore = s[-1]
    qname = s[0] 
    sname = s[1]
    q_start = s[6]
    q_end =s[7]
    s_start = min([int(i) for i in s[8:10]])
    s_end = max([int(i) for i in s[8:10]])
    qcovlen = (q_end - q_start) + 1
    scovlen = (s_end - s_start) + 1

    if len(prev) != 0:
        if qname == prev_qname and sname == prev_sname: #spanned event
            if q_start == prev_q_end:
                qcov_all += qcovlen - 1
            elif q_start > prev_q_end:
                qcov_all += qcovlen
            elif q_start <= prev_q_end:
                if q_end > prev_q_end:
                    qcov_all += q_end - prev_q_end
                else:
                    continue
            identity += blast_id
            score += bitscore
            counter += 1
            last_line = "apart"
            s_range += [s_start, s_end]
            q_range += [q_start, q_end]
        else:
            if counter == 1:
                # qcov_all = prev_q_end -  prev_q_start + 1 
                qcov_all = q_raw_cov
            mean_id = float(identity)/float(counter)
            # print s_range
            out = [prev_qname, prev_sname,mean_id, qcov_all, min(s_range), max(s_range), min(q_range), max(q_range), score ]
            out = map(lambda x:str(x), out)
            out_str = "\t".join(out)
            oot.write(out_str + "\n")
            identity = blast_id ## reset value
            score = bitscore
            counter = 1
            qcov_all = qcovlen
            s_range = [s_start, s_end]
            q_range = [q_start, q_end]

            
    else:
        s_range = [s_start, s_end]
        q_range = [q_start, q_end]
        identity = blast_id
        score = bitscore
        counter = 1
        qcov_all = qcovlen
    # print s, qcov_all, qcovlen
    #assign prevs
    prev = s
    q_raw_cov = prev[3]
    prev_s_start = s_start
    prev_s_end = s_end
    prev_q_start = q_start
    prev_q_end = q_end
    prev_id = s[2]
    prev_qname = s[0]
    prev_sname = s[1]
    prev_score = s[-1]
    

# For last line
if last_line != "apart":
    qcov_all = abs(q_start-q_end) + 1 

mean_id = float(identity)/float(counter)
out = [prev_qname, prev_sname,mean_id, qcov_all, min(s_range), max(s_range), min(q_range), max(q_range), score ]
out = map(lambda x:str(x), out)
out_str = "\t".join(out)
oot.write(out_str + "\n")
