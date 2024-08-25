'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-05-24 22:35:43
LastEditors: zpliu
LastEditTime: 2023-05-24 22:38:59
@param: 
'''
#TODO 根据与Repeat的交集结果进一步过滤SDs
from email import header
from operator import index
import pandas as pd 
import sys 

SD_filter=pd.read_csv(
    sys.argv[1],
    header=None,index_col=None,sep="\t"
)

leftSD=pd.read_csv(
    sys.argv[2],
    header=None,index_col=None,sep="\t"
)
rightSD=pd.read_csv(
    sys.argv[3],
    header=None,index_col=None,sep="\t"
)

SD_filter['leftSD']=SD_filter.apply(
    lambda x:"{}:{}-{}".format(x[0],x[1],x[2]),axis=1
)

SD_filter['rightSD']=SD_filter.apply(
    lambda x:"{}:{}-{}".format(x[3],x[4],x[5]),axis=1
)

leftSD_ID=leftSD.apply(
    lambda x:"{}:{}-{}".format(x[0],x[1],x[2]),axis=1
)
rightSD_ID=rightSD.apply(
    lambda x:"{}:{}-{}".format(x[0],x[1],x[2]),axis=1
)


#TODO 根据与TE的交集结果，筛选最终得到的SD集合
SD_filter_END=SD_filter.loc[
    (SD_filter['leftSD'].isin(rightSD_ID.values))|(SD_filter['rightSD'].isin(rightSD_ID.values))
]

#! 最终的结果
SD_filter_END[[0,1,2,3,4,5,6,7,8,9,10,11,12,13]].to_csv(sys.argv[4],header=False,index=False,sep="\t")