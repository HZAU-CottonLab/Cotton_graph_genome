'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-23 11:37:12
LastEditors: zpliu
LastEditTime: 2023-07-09 16:00:05
@param: 
'''

#TODO 根据halSynteny，提取微重排
import pandas as pd 
import sys 

def getMicroSyntenicRegion(PSLData,Syntenic_queryChr,SyntenicqueryStart,SyntenicqueryEnd,genomeType):
    #* 获取该大共线性区间内的，微共线性区域
    if genomeType=="query":
        pass 
        syntenicRegion=PSLData.loc[
             (PSLData[9]==Syntenic_queryChr)&(PSLData[11]==SyntenicqueryStart)&(PSLData[12]==SyntenicqueryEnd)
        ]
    elif genomeType=="target":
        syntenicRegion=PSLData.loc[
                (PSLData[13]==Syntenic_queryChr)&(PSLData[15]==SyntenicqueryStart)&(PSLData[16]==SyntenicqueryEnd)
        ]
    else:
        raise ("Error: please select Tag(query/target)")
    #* 以第一个共线性区间开始（得分最高）
    qStartsArray=syntenicRegion.iloc[0,19].strip(",").split(",")
    tStartsArray=syntenicRegion.iloc[0,20].strip(",").split(",")
    targetChromLength=syntenicRegion.iloc[0,14]
    AlignmentStand=syntenicRegion.iloc[0,8]
    qStartsArray=[int(i) for i in qStartsArray]
    qChrom=syntenicRegion.iloc[0,9]
    tChrom=syntenicRegion.iloc[0,13]
    qEnd=syntenicRegion.iloc[0,12]
    tEnd=syntenicRegion.iloc[0,16]
    tStart=syntenicRegion.iloc[0,15]
    Micro_syntenicRegion=[]
    if AlignmentStand=="++": 
        tStartsArray=[int(i) for i in tStartsArray]
        for index in range(0,len(tStartsArray)-1,1):
            Micro_syntenicRegion.append(
                (qChrom,qStartsArray[index],qStartsArray[index+1]-1,tChrom,tStartsArray[index],tStartsArray[index+1]-1,"++")
            )
        #* 最后一个
        Micro_syntenicRegion.append(
                (qChrom,qStartsArray[index+1],qEnd,tChrom,tStartsArray[index+1]+1,tEnd,"++")
        )
    else:
        tStartsArray=[targetChromLength-int(i)+1 for i in tStartsArray]
        # +-
        for index in range(0,len(tStartsArray)-1,1):
            Micro_syntenicRegion.append(
                (qChrom,qStartsArray[index],qStartsArray[index+1]-1,tChrom,tStartsArray[index+1]+1,tStartsArray[index],"+-")
            )
        #* 最后一个
        Micro_syntenicRegion.append(
                (qChrom,qStartsArray[index+1],qEnd,tChrom,tStart,tStartsArray[index+1],"+-")
        )
    return Micro_syntenicRegion 

if __name__=="__main__":
    PSLdata=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
    syntenicRegion=pd.read_csv(sys.argv[2],header=None,index_col=None,sep="\t")
    MicorRegionArray=pd.DataFrame()
    for Syntenic_queryChr,SyntenicqueryStart,SyntenicqueryEnd in syntenicRegion.values:
        dataArray=getMicroSyntenicRegion(PSLdata,Syntenic_queryChr,SyntenicqueryStart,SyntenicqueryEnd,"query")
        dataArray=pd.DataFrame(dataArray)
        MicorRegionArray=pd.concat([MicorRegionArray,dataArray],axis=0)
    #* 统计长度信息
    MicorRegionArray[7]=MicorRegionArray.apply(
        lambda x:x[2]-x[1]+1,axis=1
    )
    MicorRegionArray[8]=MicorRegionArray.apply(
        lambda x:x[5]-x[4]+1,axis=1
    )
    MicorRegionArray.to_csv(
        sys.argv[3],header=False,index=False,sep="\t"
    )