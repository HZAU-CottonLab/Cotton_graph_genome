'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-26 15:01:12
LastEditors: zpliu
LastEditTime: 2023-07-17 16:50:54
@param: 
'''
import sys 
import pandas as pd 
import re 
#TODO 提取断点的真实坐标：
rawVCF=pd.read_csv(sys.argv[1],header=0,index_col=None,sep="\t")

#* 每个node在对应Path中的真实坐标
#? 存在Node在对应的Path中无法对齐的情况
nodePosition=pd.read_csv(
    sys.argv[2],
    header=None,index_col=None,sep="\t"
)
#* 每个reference中的位点信息
PathName=rawVCF.columns[6:]

outData=[]
for value in rawVCF.values:
    StartNode=int(re.split(r"[><]",value[2])[1])
    EndNode=int(re.split(r"[><]",value[2])[2])
    genotypeArray=[]
    #* 省去query时间
    tmpData=nodePosition.loc[nodePosition[0].isin([StartNode,EndNode])]
    for path,genotype in zip(PathName,value[6:]):
        #* 该Path所对应的染色体编号
        Chrom=path.split("#")[1]
        startSite=tmpData.loc[(tmpData[0]==StartNode)&(tmpData[2]==path)]
        EndSite=tmpData.loc[(tmpData[0]==EndNode)&(tmpData[2]==path)]
        if startSite.empty:
            startSite=-1
        else:
            startSite=startSite.iloc[0,1]
        if EndSite.empty:
            EndSite=-1 
        else:
            EndSite=EndSite.iloc[0,1]
        #* 每个材料中的真实位置
        genotypeArray.append(
            "{},{}:{}-{}".format(genotype,Chrom,startSite,EndSite)
        )
    #* SVs在参考基因组中对应的坐标
    refGenotype,refChr,refSite1,refSite2=re.findall(r"([^,]+),([^:]+):(-?[^-]+)-(-?[0-9NA]+)",genotypeArray[0])[0]
    if refSite1=="-1" or refSite2=="-1":
        #* SVs断点在参考基因组中的坐标
        refstart=value[1]
        refEnd=refstart+value[3]
    else:
        #* VG paths得到
        refstart=min(int(refSite1),int(refSite2))
        refEnd=max(int(refSite1),int(refSite2))

    outData.append(
        [refChr,refstart,refEnd]+list(value[2:6])+genotypeArray
    )
#* 填补好SVs断点真实坐标的最终文件
outData=pd.DataFrame(outData,columns=['#Chr','start','end']+list(rawVCF.columns)[2:])

outData.to_csv(sys.argv[3],header=True,index=False,sep="\t")