'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-28 17:50:39
LastEditors: zpliu
LastEditTime: 2023-07-08 21:38:11
@param: 
'''
from email import header
import pandas as pd 
import sys 
filter_PS=pd.read_csv(
    sys.argv[1],header=0,index_col=None,sep="\t"
)
rawPAV=pd.read_csv(
    sys.argv[2],header=None,index_col=None,sep="\t"
)

filter_PS.index=filter_PS.apply(
    lambda x:"{}_{}".format(x['#Chr'],x['ID']),axis=1
)

rawPAV['SV_id']=rawPAV.apply(
    lambda x:"{}_{}".format(x[1],x[0]),axis=1
)



#TODO  不同的染色体可能有相同的SVid，
filterData=rawPAV.loc[rawPAV['SV_id'].isin(
    filter_PS.index
)]

filterData=filterData[[1,2,3,4,5,6,7,8]]
#* 过滤掉-1项目
filterData=filterData.loc[
    (filterData[2]!=-1)&(
        filterData[3]!=-1
    )&(
        filterData[5]!=-1
    )&(
        filterData[6]!=-1
    )
]
filterData.to_csv(sys.argv[3],header=False,index=False,sep="\t")