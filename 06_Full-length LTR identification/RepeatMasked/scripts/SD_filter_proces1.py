'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-24 22:35:32
LastEditors: zpliu
LastEditTime: 2023-05-24 22:36:33
@param: 
'''
#TODO 对原始的SD序列根据长度、相似度以及gap率进行过滤
import pandas as pd 
import sys
SD_annotion=pd.read_csv(
   sys.argv[1],
    header=None,index_col=None,sep="\t"
)

def SD_filter(SD):
    X,ID=SD[13].split(";")
    Xratio=float(X.strip("X="))
    IDratio=float(ID.strip("ID="))
    if SD[11]>=1000 and Xratio<=10 and IDratio<=50:
        #* 片段长度超过1000
        #* 相似度超过90%
        #* gap少于50%
        return True
    else:
        return False


SD_filter_Index=SD_annotion.apply(lambda x:SD_filter(x),axis=1)
SD_filter=SD_annotion.loc[SD_filter_Index]


SD_filter.to_csv(sys.argv[2],header=False,index=False,sep="\t")