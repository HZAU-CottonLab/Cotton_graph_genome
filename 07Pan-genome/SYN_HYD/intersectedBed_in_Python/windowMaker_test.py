'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-07-15 23:04:47
LastEditors: zpliu
LastEditTime: 2023-07-16 15:37:04
@param: 
'''
from bedtools import windowMaker 
import pandas as pd  

A2_geneBed=pd.read_csv(
    "J85_gene.bed",
    header=None,index_col=None,sep="\t"
)


#TODO 三种窗口
geneBody=[]
TSS=[]
TTS=[]
for val in A2_geneBed.values:
    geneBody.append(val)
    if val[4]=="+":
        TTS.append(
                (val[0],val[2]+1,val[2]+2000,val[3],val[4])
        )
        if val[1]>2000:
            TSS.append(
                (val[0],val[1]-2000,val[1]-1,val[3],val[4])
            )
        else:
            TSS.append(
                (val[0],1,val[1]-1,val[3],val[4])
            )
    else:
        #* 负链
        TSS.append(
                (val[0],val[2]+1,val[2]+2000,val[3],val[4])
        )
        if val[1]>2000:
            TTS.append(
                (val[0],val[1]-2000,val[1]-1,val[3],val[4])
            ) 
        else:
            TTS.append(
                (val[0],1,val[1]-1,val[3],val[4])
            )
#* 三种区间制作windows
TSS=pd.DataFrame(TSS)
geneBody=pd.DataFrame(geneBody)
TTS=pd.DataFrame(TTS)


#* 对不同的链的坐标进行处理
TTS_window=pd.DataFrame()
for stand,standItem in TTS.groupby([4]):
    if stand=="+":
        splitData=windowMaker(standItem,seqRevered=True,seqAnnotation="TTS",windowSize=200,stepSize=100)
        TTS_window=pd.concat([TTS_window,splitData],axis=0)
    else:
        splitData=windowMaker(standItem,seqRevered=False,seqAnnotation="TTS",windowSize=200,stepSize=100)
        TTS_window=pd.concat([TTS_window,splitData],axis=0) 

TSS_window=pd.DataFrame()
for stand,standItem in TSS.groupby([4]):
    if stand=="+":
        splitData=windowMaker(standItem,seqRevered=True,seqAnnotation="TSS",windowSize=200,stepSize=100)
        TSS_window=pd.concat([TSS_window,splitData],axis=0)
    else:
        splitData=windowMaker(standItem,seqRevered=False,seqAnnotation="TSS",windowSize=200,stepSize=100)
        TSS_window=pd.concat([TSS_window,splitData],axis=0) 

gene_window=pd.DataFrame()
for stand,standItem in geneBody.groupby([4]):
    if stand=="+":
        splitData=windowMaker(standItem,seqRevered=True,seqAnnotation="gene",windowNum=200)
        gene_window=pd.concat([gene_window,splitData],axis=0)
    else:
        splitData=windowMaker(standItem,seqRevered=False,seqAnnotation="gene",windowNum=200)
        gene_window=pd.concat([gene_window,splitData],axis=0) 

#* 同一个基因ID，合并所有的windows ID, 并使的windows ID连续

