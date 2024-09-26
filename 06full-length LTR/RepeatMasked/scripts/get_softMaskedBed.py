'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-22 10:21:15
LastEditors: zpliu
LastEditTime: 2023-05-22 10:21:16
@param: 
'''
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO.FastaIO import FastaWriter
import pandas as pd
import re
import sys 
#TODO 提取soft-masked基因组中的Repeat的坐标
RepeatIndex = []
with open(sys.argv[1]) as handle:
    for record in FastaIterator(handle):
        print("deal with {}".format(record.id))
        for f in re.finditer('[actg]+',str(record.seq)):
            start,end=f.span()
            RepeatIndex.append(
                (
                    record.id,start+1,end
                )
            )

pd.DataFrame(RepeatIndex).to_csv(sys.argv[2],header=False,index=False,sep="\t")
                