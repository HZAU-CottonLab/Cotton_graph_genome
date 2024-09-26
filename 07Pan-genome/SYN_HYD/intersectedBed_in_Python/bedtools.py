'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-07-11 19:37:45
LastEditors: zpliu
LastEditTime: 2023-07-18 16:37:28
@param: 
'''
import pandas as pd 
import pybedtools 
from tempfile import NamedTemporaryFile


#TODO 根据基因组共线性区间，提取高度差异的区域
def intersectBed(queryDataFram,targetDataFram):
    #* 将两个bed文件取交集
    queryBedFile=NamedTemporaryFile(mode='w+')
    targetBedFile=NamedTemporaryFile(mode='w+')
    queryDataFram.to_csv(
        queryBedFile.name,header=False,index=False,sep="\t"
    )
    targetDataFram.to_csv(
        targetBedFile.name,header=False,index=False,sep="\t"
    )
    queryBed=pybedtools.BedTool(queryBedFile.name)
    targetBed=pybedtools.BedTool(targetBedFile.name)
    intersectBed=queryBed.intersect(
        targetBed,loj=True
    )
    out=[]
    for val in intersectBed:
        out.append(
            str(val).strip("\n").split("\t")
        )
    out=pd.DataFrame(out)
    # out=out.astype(
    #     {
    #         1:int,
    #         2:int,
    #         4:int,
    #         5:int
    #     }
    # )
    queryBedFile.close()
    targetBedFile.close()
    return out  


def merge_bed(DataFrame):
    '''
    #* 合并基因组坐标
    '''
    bedFile=NamedTemporaryFile(mode="w+")
    DataFrame.to_csv(bedFile.name,header=False,index=False,sep="\t")
    input_bed=pybedtools.BedTool(bedFile)
    mergeBed=input_bed.merge(d=0)
    out=[]
    for i in mergeBed:
        out.append(
            str(i).strip("\n").split("\t")
        )
    out=pd.DataFrame(out)
    out=out.astype(
        {
            1:int,
            2:int
        }
    )
    bedFile.close()
    return out


def subtractBed_genome(fastaFaiFile,inputBedDataFrame):
    fastaFaiDataFrame=pd.read_csv(
        fastaFaiFile,header=None,index_col=None,sep="\t"
    )
    #* 基因组文件
    fastaBedFile=NamedTemporaryFile(mode="w+")
    fastaFaiDataFrame[5]=1
    fastaFaiDataFrame[[0,5,1]].to_csv(
        fastaBedFile.name,header=False,index=False,sep="\t"
    )
    #* 需要减去的坐标
    inputBedFile=NamedTemporaryFile(mode="w+")
    inputBedDataFrame.to_csv(inputBedFile.name,header=False,index=False,sep="\t")
    fastaBed=pybedtools.BedTool(fastaBedFile.name)
    inputBed=pybedtools.BedTool(inputBedFile.name)
    substractBed=fastaBed.subtract(
        inputBed
    )
    out=[]
    # print(substractBed.head())
    for val in substractBed:
        out.append(
            str(val).strip("\n").split("\t")
        )
    out=pd.DataFrame(out)
    out=out.astype(
        {
            1:int,
            2:int
        }
    )
    inputBedFile.close()
    fastaBedFile.close()
    return out 


def windowMaker(dataFrame,seqRevered,**kwargs):
    '''
    @dataFrame: BED
    @seqRevered: sequence reverse or not
    @**kwargs
    @windowNum: count of window 
    @seqAnnotation: seq location, gene, upstream
    @windowSize: fixed windowSize
    @stepSize: sliding windowSize
    '''
    #* 将bed拆分为多个windows
    inputBedFile=NamedTemporaryFile(mode="w+")
    dataFrame.to_csv(inputBedFile.name,header=False,index=False,sep="\t")
    inputBed=pybedtools.BedTool(inputBedFile.name)
    if kwargs.get("windowNum"):
        windowSplit=inputBed.window_maker(
            b=inputBedFile.name,
            reverse=seqRevered,
            i="srcwinnum",
            n=kwargs.get("windowNum"),
        )
    else:
        windowSplit=inputBed.window_maker(
            b=inputBedFile.name,
            reverse=seqRevered,
            i="srcwinnum",
            w=kwargs.get("windowSize"),
            s=kwargs.get("stepSize"),
        ) 
    seqAnnotation=kwargs.get("seqAnnotation")
    out=[]
    for val in windowSplit:
        val=str(val).strip("\n").split("\t")
        #* 添加ID和window ID以及注释信息
        out.append(
            (
                val[0],int(val[1]),
                int(val[2]),"_".join(val[3].split("_")[0:-1]),
                int(val[3].split("_")[-1]),seqAnnotation
            )
        )
    inputBedFile.close()
    return pd.DataFrame(out)