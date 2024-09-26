'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-06-26 15:01:12
LastEditors: zpliu
LastEditTime: 2024-09-26 09:42:25
@param: 
'''
import sys 
import pandas as pd 
import re 
rawVCF=pd.read_csv(sys.argv[1],header=0,index_col=None,sep="\t")

#* The real coordinates of each node in the corresponding Path.
#? There is a situation that the Node cannot be aligned in the corresponding Path.
nodePosition=pd.read_csv(
    sys.argv[2],
    header=None,index_col=None,sep="\t"
)
#* Site information in each reference
PathName=rawVCF.columns[6:]

outData=[]
for value in rawVCF.values:
    StartNode=int(re.split(r"[><]",value[2])[1])
    EndNode=int(re.split(r"[><]",value[2])[2])
    genotypeArray=[]
    tmpData=nodePosition.loc[nodePosition[0].isin([StartNode,EndNode])]
    for path,genotype in zip(PathName,value[6:]):
        #* Chromosome number corresponding to this Path.
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
        #* The true position in each material
        genotypeArray.append(
            "{},{}:{}-{}".format(genotype,Chrom,startSite,EndSite)
        )
    #* Coordinates of SVs in the reference genome
    refGenotype,refChr,refSite1,refSite2=re.findall(r"([^,]+),([^:]+):(-?[^-]+)-(-?[0-9NA]+)",genotypeArray[0])[0]
    if refSite1=="-1" or refSite2=="-1":
        #* Coordinates of SVs breakpoint in reference genome
        refstart=value[1]
        refEnd=refstart+value[3]
    else:
        #* VG paths
        refstart=min(int(refSite1),int(refSite2))
        refEnd=max(int(refSite1),int(refSite2))

    outData.append(
        [refChr,refstart,refEnd]+list(value[2:6])+genotypeArray
    )
#* Fill in the final file of the real coordinates of SVs breakpoint.
outData=pd.DataFrame(outData,columns=['#Chr','start','end']+list(rawVCF.columns)[2:])

outData.to_csv(sys.argv[3],header=True,index=False,sep="\t")