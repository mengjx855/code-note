#!/share/data2/guorc/Software/conda/py3/bin/python
# -*- coding: UTF-8 -*-

import os,sys,argparse
def code_help():
    
    args = sys.argv
    script_path = os.path.abspath(args[0])

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--ffn", type=str, help="ffn path from Prodigal")
    parser.add_argument("-b", "--vcf", type=str, help="vcf file from bcftools")
    parser.add_argument("-o", "--outdir", type=str, help="name of output dir")
    parser.add_argument("-d", "--minDepth", type=int, help="minimum depth, default 40", default=40)
    parser.add_argument("-c", "--SNPcoverage", type=float,help="SNP coverage, default 0.05", default=0.05)
    parser.add_argument("-t", "--threads", type=int, help="threads, default 10", default=10)
    parser.add_argument("-s", "--snp", action='store_true', help="only count the number of SNPs", default=False)
    args = parser.parse_args()

    if not args.ffn or not args.vcf or not args.outdir:
        os.system('python %s -h'%script_path)
        sys.exit()
    
    #if os.path.isfile(args.outfile):
    #    print('\nExit: ' + args.outfile + ' is exist!!\n')
    #    sys.exit()
    
    return args.ffn,args.vcf,args.outdir,args.minDepth,args.SNPcoverage,args.threads,args.snp

ffnPath, vcfPath, outdir, minDepth, SNPcoverage, Threads, SNP = code_help()

#================================================================================================

import re
import pandas as pd
from math import ceil
import numpy as np
from multiprocessing import Process
import time
import gc
import itertools

# ffnPath='00.data/ref.ffn'
# vcfPath='03.snp_grad/KD-243.vcf.filter'
# outdir='aa'
# minDepth=40
# SNPcoverage=0.05
# Threads = 10
# SNP=False

#Data Preparation
#codon11Path='/share/data5/guorc/database/nt_20190330/codon.11_20191031'
codonNSPath='/share/data2/guorc/script/WGS/Variation/pNpS_data/sn_tot.list1'
#codonNStranPath='/share/data2/guorc/script/WGS/Variation/pNpS_data/sn_tot.list2'


def openf(Path,Line=False):
    with open(Path) as f:
        if Line:
            inputFile=f.readlines()
        else:
            inputFile=f.read()
    return inputFile

def ffnTreat(fPath):
    inputFile = openf(fPath)
    ffn = inputFile.split('>')[1:]
    dic = {}
    for x in ffn:
        Tmp = x.split('\n',1)
        Seq = Tmp[1].replace('\n','')
        Tmp = Tmp[0].split(' ')
        ContigID = Tmp[0].rsplit('_',1)[0]
        if ContigID in dic:
            dic[ContigID].append({'SeqID':Tmp[0],'Start':int(Tmp[2]),\
            'End':int(Tmp[4]),'Direction':Tmp[6],'Seq':Seq})
        else:
            dic[ContigID] = [{'SeqID':Tmp[0],'Start':int(Tmp[2]),\
            'End':int(Tmp[4]),'Direction':Tmp[6],'Seq':Seq}]
    return dic

def vcfTreat(vPath,minD=False,SNPcov=False):
    vcf = vPath#[vPath[0].isin(list(ffn.keys()))]
    I16 = [list(map(int,list(re.search(':(\d+),(\d+),(\d+),(\d+)$',x).groups()))) for x in vcf[9]]
    vcf[10] = [sum(x) for x in I16]
    Filter = [sum(x)>=minD for x in I16]
    vcf = vcf[Filter]
    I16 = [x for x in I16 if sum(x)>=minD]
    cumulative_depth = sum(list(itertools.chain.from_iterable(I16)))
    Filter=[]
    for x in I16:
        if sum(x[2:])<=sum(x[:2]):
            Filter.append(False)
        else:
            Filter.append(True)
    vcf.loc[Filter,3],vcf.loc[Filter,4]=vcf.loc[Filter,4],vcf.loc[Filter,3]
    Filter=[]
    for x in I16:
        if sum(x[2:])/sum(x)>=SNPcov and sum(x[2:])/sum(x)<(1-SNPcov):
            Filter.append(False)
        else:
            Filter.append(True)
    vcf.loc[Filter,4]=vcf.loc[Filter,3]
    vcf.loc[vcf[3]==vcf[4],10] = np.nan
    return vcf,cumulative_depth

def Complement(Base):
    Base_C = {'A':'T','T':'A','G':'C','C':'G'}
    Base = Base_C[Base]
    return Base

def getSNPcodon(vcf,ffn):
    vcf_f=vcf[[0,1,3,4,10]].reset_index(drop=True)
    vcf_f.columns = ['ContigID','POS','Ref','Aln','SNPposDepth']
    vcf_f = vcf_f[vcf_f['ContigID'].isin(list(ffn.keys()))]
    Refcodon=[]
    SNPcodon=[]
    SeqID=[]
    Ncodon=[]
    Spos=[]
    Epos=[]
    SNPcodon = pd.DataFrame()
    for x,y in ffn.items():
        tab = vcf_f.loc[vcf_f['ContigID']==x]
        if tab.shape[0]==0:
            continue
        for z in y:
            tabf = tab[(tab['POS']>=z['Start']) & (tab['POS']<=z['End'])].reset_index(drop=True)
            if tabf.shape[0]==0:
                continue
            
            #Ref_codon = re.findall(r'.{3}',z['Seq'])
            SNP_seq = list(z['Seq'])
            Ref_seq = list(z['Seq'])
            SNPposDepth = [np.nan for x in range(len(z['Seq']))]
            POS = [i-z['Start']+1 for i in tabf['POS']]
            Ncodon = [ceil(j/3) for j in POS]
            
            if z['Direction']=='1':
                Spos = [z['Start']+3*(j-1) for j in Ncodon]
                Epos = [j+2 for j in Spos]
                for i in range(len(POS)):
                    SNP_seq[POS[i]-1] = tabf.loc[i,'Aln'][0]
                    Ref_seq[POS[i]-1] = tabf.loc[i,'Ref'][0]
                    SNPposDepth[POS[i]-1] = tabf.loc[i,'SNPposDepth']
                SNP_codon = re.findall(r'.{3}',''.join(SNP_seq))
                SNP_codon = [SNP_codon[x-1] for x in Ncodon]
                Ref_codon = re.findall(r'.{3}',''.join(Ref_seq))
                Ref_codon = [Ref_codon[x-1] for x in Ncodon]
                #SNPposDepth = re.findall(r'.{3}',''.join(SNPposDepth))
                SNPposDepth = [[j for j in SNPposDepth[i:i+3] if not np.isnan(j)] for i in range(0, len(SNPposDepth), 3)]
                SNPposDepth = [','.join(map(str,SNPposDepth[x-1])) for x in Ncodon]
            else:
                Epos = [z['Start']+3*(j-1) for j in Ncodon]
                Spos = [j+2 for j in Epos]
                for i in range(len(POS)):
                    Base=Complement(tabf.iloc[i,3][0]) #3==Aln
                    SNP_seq[-POS[i]] = Base
                    Base=Complement(tabf.iloc[i,2][0])
                    Ref_seq[-POS[i]] = Base
                    SNPposDepth[-POS[i]] = tabf.iloc[i,4]
                SNP_codon = re.findall(r'.{3}',''.join(SNP_seq))
                SNP_codon = [SNP_codon[-x] for x in Ncodon]
                Ref_codon = re.findall(r'.{3}',''.join(Ref_seq))
                Ref_codon = [Ref_codon[-x] for x in Ncodon]
                SNPposDepth = [[j for j in SNPposDepth[i:i+3] if not np.isnan(j)] for i in range(0, len(SNPposDepth), 3)]
                SNPposDepth = [','.join(map(str,SNPposDepth[-x])) for x in Ncodon]
                Ncodon = [int(len(z['Seq'])/3-(x-1)) for x in Ncodon]
            
            Tab = pd.DataFrame({'#SeqID':z['SeqID'],'Refcodon':Ref_codon,'SNPcodon':SNP_codon,'Ncodon':Ncodon,'Start':Spos,'End':Epos,'SNPposDepth':SNPposDepth})
            Filter = Tab['Ncodon'].value_counts()
            Tab = Tab[Tab['Ncodon'].isin(Filter.index[Filter.values==3])].drop_duplicates()
            Tab['Ncodon'] = Tab['Ncodon']-1
            SNPcodon = pd.concat([SNPcodon, Tab], axis=0)
    return SNPcodon

def SNPCombType(Ref,Aln):
    '''
    发生突变的codon，生成独立的突变codon组合
    '''
    A = [Aln[0]+Ref[1]+Ref[2], Ref[0]+Aln[1]+Ref[2], Ref[0]+Ref[1]+Aln[2]]
    return A

def getpNpS(SNPcodon,bPath):
    codon_NS = pd.read_csv(bPath,sep='\t',header=0)
    codon_NS.columns = ['Refcodon','N','S']
    codon_NS['AA'] = [x[5] for x in codon_NS['Refcodon']]
    codon_NS['Refcodon'] = [x[0:3] for x in codon_NS['Refcodon']]
    SNPcodon = pd.merge(SNPcodon,codon_NS,how='left',on='Refcodon')
    #统计同义和非同义
    SNPcodon['Nd'] = 0; SNPcodon['Sd'] = 0
    Poly = SNPcodon[SNPcodon['Refcodon']!=SNPcodon['SNPcodon']].reset_index(drop=True)
    for x in range(Poly.shape[0]):
        SNP = [i for i in list(set(SNPCombType(Poly.iloc[x,1],Poly.iloc[x,2]))) if i!=Poly.iloc[x,1]]
        for y in range(len(SNP)):
            AA = codon_NS['AA'][codon_NS['Refcodon'] ==SNP[y]].values[0]
            if AA == Poly.loc[x,'AA']:
               Poly.loc[x,'Sd'] = Poly.loc[x,'Sd'] + 1
            else:
               Poly.loc[x,'Nd'] = Poly.loc[x,'Nd'] + 1
    SNPcodon = pd.concat([Poly,SNPcodon[SNPcodon['Refcodon']==SNPcodon['SNPcodon']]], axis=0)
    SNPcodon.drop(['AA'],axis=1,inplace=True)
    #基因组的pN_pS计算
    #pN_pS = pN_pS_1(SNPcodon['Nd'],SNPcodon['N'],SNPcodon['Sd'], SNPcodon['S'])
    pN_pS = pN_pS_2(SNPcodon['Nd'],SNPcodon['N'],SNPcodon['Sd'], SNPcodon['S'], SNPcodon['SNPposDepth'])
    #基因的pN_pS计算
    GeneData = SNPcodon.groupby(['#SeqID'])
    Genecodon = {}
    for x in list(GeneData.groups):
        Gene_data = GeneData.get_group(x)
        if Gene_data.shape[0] >= 3: #满足测序深度的总codon太少对于一个基因而言，计算的pNpS不具备参考意义
            Genecodon[x]=pN_pS_2(Gene_data['Nd'],Gene_data['N'],Gene_data['Sd'], Gene_data['S'], Gene_data['SNPposDepth'])
    Genecodon = [key+"\t"+str(value)+"\n" for key, value in Genecodon.items() if value != 'NA']
    return pN_pS,SNPcodon,Genecodon

def pN_pS_1(Nd,N,Sd,S):
    #pN_pS计算，添加一个极小值，避开分母为0
    pseudo_count = 1e-20
    if sum(Nd)==0 and sum(Sd)==0:
        pN_pS = 'NA'
    else:
        pN_pS = ((sum(Nd)+pseudo_count)/(sum(N)+pseudo_count))/((sum(Sd)+pseudo_count)/(sum(S)+pseudo_count))
    return pN_pS

def pN_pS_2(Nd,N,Sd,S,SNPdepth):
    #参考文献：https://peerj.com/articles/2959
    import statistics, math
    if sum(Nd)==0 and sum(Sd)==0:
        pN_pS = 'NA'
    else:
        #pN_pS计算，添加snp位点测序深度中位数的开方的一半，校正测序深度和避开分母为0 
        pseudo_count = [float(item.strip()) for sublist in SNPdepth[SNPdepth!=''].str.split(',') for item in sublist]
        pseudo_count = math.sqrt(statistics.median(pseudo_count))/2
        pN_pS = ((sum(Nd)+pseudo_count)/sum(N))/((sum(Sd)+pseudo_count)/sum(S))
    return pN_pS


def run(buk,Group,Ffn,minDepth,SNPcoverage,bPath,SNP):
    for group in buk:
        #group=N[9][3]
        #time.ctime(time.time())
        vcf = Group.get_group(group).drop([0],axis=1).reset_index(drop=True)
        vcf.columns = range(vcf.shape[1])
        vcf[1] = vcf[1].astype("int")
        ff = vcf[0].drop_duplicates().tolist()
        ffn = dict((key,value) for key,value in Ffn.items() if key in ff)
        vcf, cumulative_depth = vcfTreat(vcf,minDepth,SNPcoverage)
        
        #计算snp个数及满足测序深度的位点数及这些位点的累计测序深度
        with open(outdir+'/'+group+'.snp','w') as f:
            f.write('#SNP Sites: '+str(vcf[vcf.iloc[:,3]!=vcf.iloc[:,4]].shape[0])+' '+str(vcf.shape[0])+' '+str(cumulative_depth)+'\n')
        vcf.loc[vcf.iloc[:,3]!=vcf.iloc[:,4],0:1].to_csv(outdir+'/'+group+'.snp',sep='\t',index=0, header=0, mode='a')
        if SNP:
            continue
        
        #####
        SNPcodon = getSNPcodon(vcf,ffn)
        if SNPcodon.shape[0] < 3: #满足测序深度的总codon太少对于一个基因组而言，计算的pNpS不具备参考意义
            print('\n'+outdir+'/'+group+': The number of detectable bases is too litte!!!!\n')
            continue
        
        pN_pS,SNPcodon,Genecodon = getpNpS(SNPcodon,bPath)
        with open(outdir+'/'+group+'.pnps','w') as f:
            f.write('#pN/pS = '+str(pN_pS)+'\n')
            f.writelines(Genecodon)
        SNPcodon.to_csv(outdir+'/'+group+'.pnps_data',sep='\t',index=False)

##################################################################
#absPath=os.path.abspath(ffnPath)
vPath=pd.read_csv(vcfPath,sep='\t',header=None,dtype=object)
Group = vPath.groupby([0])
del vPath
gc.collect()

Ffn = ffnTreat(ffnPath)
N = list(Group.groups)
N = [N[x:x+Threads] for x in range(0,len(N),Threads)]

threads = []
for x in N:
    t = Process(target=run, args=(x,Group,Ffn,minDepth,SNPcoverage,codonNSPath,SNP))
    threads.append(t)
    t.start()

for x in threads:
    x.join()

