'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-28 17:26:53
LastEditors: zpliu
LastEditTime: 2023-06-28 17:26:54
@param: 
'''
from email import header
import sys 
import pandas as pd 
#TODO 过滤parent snarl
rawPAV=pd.read_csv(
    sys.argv[1],header=0,index_col=None,sep="\t"
)


parent_snarl=[]
filterData=[]
for index in range(0,rawPAV.shape[0]-1,1):
    nextChr,nexStart,nextEnd=rawPAV.iloc[index+1,0:3]
    SVChr,SVstart,SVend=rawPAV.iloc[index,0:3]
    if SVChr==nextChr and SVstart<=nexStart and SVend>=nextEnd:
        #* 该SVs是parent snarl
        parent_snarl.append(
            rawPAV.iloc[index,:]
        ) 
    else:
        filterData.append(
            rawPAV.iloc[index,:]
        ) 
filterData=pd.DataFrame(filterData)
parent_snarl=pd.DataFrame(parent_snarl)

filterData.to_csv(sys.argv[2],header=True,index=False,sep="\t")
parent_snarl.to_csv(sys.argv[3],header=True,index=False,sep="\t")