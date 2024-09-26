'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-09-01 19:57:25
LastEditors: zpliu
LastEditTime: 2023-10-20 21:52:58
@param: 
'''
import pandas as pd 
import sys 
import re 
import logging
logging.basicConfig(
    level=logging.INFO
)
#* 导入所需要的包
#* 导入所需要的包
sys.path.append("/public/home/zpliu/Pan-genome/SV_parallele_V2/GenomeCompare")
from divergenceRegion.bedtools import intersectBed 


#! 和syntenic区域同样的操作
A2_sampleInfo = pd.read_csv(
    "/public/home/zpliu/Pan-genome/SV_parallele_V2/superPan/A2_longRead_call/15.Assembly.txt",
    header=None, index_col=None, sep="\t"
)
AD1_sampleInfo = pd.read_csv(
    "/public/home/zpliu/Pan-genome/SV_parallele_V2/wild_cultivar/AD1_35_info.txt",
    header=None, index_col=None, sep="\t"
)
geneOrthoBedFile={}
AllSample=pd.concat(
    [A2_sampleInfo[0],AD1_sampleInfo[0]],axis=0
)
# for sample in AD1_sampleInfo[0].values:
#     geneOrthoBedFile[sample]=geneOrthoBedFile.get(sample, 
#         pd.read_csv(
#                 "/public/home/zpliu/Pan-genome/SV_parallele_V2/convergnce_divergence/JCVI/A2_At_Dt_D5/{}.Dt_gene_orthogroup.txt".format(
#                     sample
#                 ),
#                 header=None,index_col=None,sep="\t"
#             )
#     )
#? D5材料的同源基因    
geneOrthoBedFile["D5"]=geneOrthoBedFile.get("D5", 
        pd.read_csv(
                "/public/home/zpliu/Pan-genome/SV_parallele_V2/convergnce_divergence/JCVI/A2_At_Dt_D5/{}.pep_gene_orthogroup.txt".format(
                    "D5"
                ),
                header=None,index_col=None,sep="\t"
            )
    )

#* 获取syntenic区间内，同源基因在所有材料间的对应关系
def get_queryRegion_OrthoId(queryRegion:pd.DataFrame,geneBed:pd.DataFrame,genome:str):
    '''
    @description 统计指定区间内,同源基因的ID信息
    '''
    intersectedData=intersectBed(
            queryRegion,geneBed
        )
    #!获得完全位于区域内的基因IDs
    filter_intersected=intersectedData.loc[ 
            intersectedData.apply(lambda x:True if int(x[4])>=int(x[1]) and int(x[5])<=int(x[2]) else False,axis=1)
        ]
    if filter_intersected.empty:
        #* 该区域内没有基因信息
        return pd.DataFrame(columns=["orthoId",genome]) 
    #! 对于单个区域内同时存在多个相同的Orthogroup,需要对这个orthogroup进行编号
    reOrderData=[]
    for OrthoGroupId,GroupData in filter_intersected.groupby([8]):
        reOrderData.append(
            (OrthoGroupId,"-".join(GroupData[6].values))
        )
    reOrderData=pd.DataFrame(reOrderData,columns=["orthoId",genome])
    return  reOrderData


#! 统计每个syntenic内同源基因在50个基因组中的情况
syntenicGeneInfo=[]
sampleOrder=[ 'D5']

At_divergenceRegion=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
for val in At_divergenceRegion.values:
    #* 获取该区域内所有材料的同源基因对应情况
    logging.info(
        "deal with: >>>>>>>>>>>>>>>>>>>\n{}".format(
            "\t".join([str(i) for i in val])
        )
    )
    regionGenome=pd.DataFrame(columns=["orthoId"])
    for sample,genomeRegion in zip(sampleOrder,val[1:]):
        if genomeRegion==".":
            sampleRegion=pd.DataFrame(
                columns=["orthoId",sample]
            )
        else:
            #*非空的基因组区间
            Chr=genomeRegion.split(":")[0]
            Start=int(genomeRegion.split(":")[1].split("-")[0])
            End=int(genomeRegion.split("-")[-1])
            #! 与基因Orthogroup取交集
            geneBedFile=geneOrthoBedFile.get(
                sample
            )
            sampleRegion=get_queryRegion_OrthoId(
                pd.DataFrame(
                    [(Chr,Start,End)]
                ),
                geneBedFile,
                sample
            )
        #* 合并该区域在所有材料的同源基因信息
        regionGenome=pd.merge(
            regionGenome,
            sampleRegion,
            on=['orthoId'],
            how='outer'
        )
    #! 调整merge后的顺序
    regionGenome=regionGenome[['orthoId']+sampleOrder]
    for geneInfo in  regionGenome.values:
        syntenicGeneInfo.append(
            list(val)+ list(geneInfo)
        )
    
syntenicGeneInfo=pd.DataFrame(syntenicGeneInfo,
   columns=   [
         'windowsType','D5','orthId','D5_gene'
   ]
 )
#! 有些divergence区间在35个材料间都不存在基因
syntenicGeneInfo.to_csv(
    sys.argv[2],header=True,index=False,sep="\t",
    na_rep="."
)