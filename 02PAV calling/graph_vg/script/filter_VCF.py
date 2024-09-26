'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-25 19:47:12
LastEditors: zpliu
LastEditTime: 2024-09-26 09:44:48
@param: 
'''
import sys 
import pandas as pd 
#TODO Extract PAV information from vcf file obtained by vg deconstruct.
#* Get the ID of breakpoint in graph, and get its coordinates in genome according to vg find.
skipRows=int(sys.argv[1])
rawVCF=pd.read_csv(
    sys.argv[2],
    header=0,index_col=None,sep="\t",skiprows=skipRows
)
referencePath=sys.argv[3]

filterValue=[]
for value in rawVCF.values:
    #TODO 提取长度差超过50bp的PAVs
    refSeqLength=len(value[3])
    AltSeqLengthArray=[len(val) for val in value[4].split(",")]
    #* The maximum length is more than 50bp.
    maxSeqDivergence=max([abs(refSeqLength-i) for i in AltSeqLengthArray])
    #* SV levels(0=top level)
    SV_level=value[7].split(";")[-1]
    if maxSeqDivergence>=50:
        SV_info=[value[0],value[1],value[2],refSeqLength,",".join([str(i) for i in AltSeqLengthArray]),SV_level]
        genotype=[0]+list(value[9:])
        mergeInfo=SV_info+genotype
        filterValue.append( 
           mergeInfo
        ) 
    else:
        #* 较小的结构变异?
        continue  
filterValue=pd.DataFrame( 
    filterValue,columns=["#Chr","Start",'ID','refLen','AltLen',"SV_level"]+[referencePath]+list(rawVCF.columns[9:])
)

filterValue.to_csv(
    sys.argv[4],header=True,index=False,sep="\t"
)