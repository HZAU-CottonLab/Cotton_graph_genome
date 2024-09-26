'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-10-21 21:10:31
LastEditors: zpliu
LastEditTime: 2023-10-23 16:02:19
@param: 
#! 基于四个基因组中模糊的共线性区间,获取对应的pan-gene
'''
import pandas as pd 
import sys 
import logging
logging.basicConfig(
    level=logging.INFO
)

#* 每个基因组在所有材料的OrthoId信息
logging.info(
    "reading liftover gene File..."
)

At_synGene=pd.read_csv(
    "At_all_window_liftover_35Gene.txt",
    header=0,index_col=None,sep="\t"
)
Dt_synGene=pd.read_csv(
    "Dt_all_window_liftover_35Gene.txt",
    header=0,index_col=None,sep="\t"
)
A2_synGene=pd.read_csv(
    "A2_all_window_liftover_15Gene.txt",
    header=0,index_col=None,sep="\t"
)
#* D5单个材料的gene统计
D5_synGene=pd.read_csv(
    "D5_all_window_Gene.txt",
    header=0,index_col=None,sep="\t"
)

# A2_D5_queryWindows=pd.read_csv(
#     "A2_D5_mosaic_SYN.txt",
#     header=None,index_col=None,sep="\t"
# )
A2_D5_queryWindows=pd.read_csv(
    sys.argv[1],
    header=None,index_col=None,sep="\t"
)
def mergeOrthoGroupData(geneData,newColumns):
    '''
    由于mosica区间内存在多个windows,因此合并多个windows中的同源基因Group
    '''
    if geneData.empty:
        #! mosaic区域在该区域没有对应的交集
        return pd.DataFrame(columns=newColumns)
    outData=[]
    for OrthoId,OrthoIdData in geneData.groupby(['orthId']):
            mergeOrthGroup=[ 
                OrthoId
            ]
            for sample,sampledata in OrthoIdData.iloc[:,1:].iteritems():
                OrthoIDGene=[i for i in sampledata if i!="."] 
                if len(OrthoIDGene)==0:
                    mergeOrthGroup.append(".")
                else:
                    mergeOrthGroup.append("-".join(OrthoIDGene))
            outData.append(
                mergeOrthGroup
            )
    outData=pd.DataFrame(outData,columns=newColumns)
    return outData



GeneData=[]
for value in A2_D5_queryWindows.values:
    logging.info(
        "deal with >>>>>>\n{}:{}-{} ".format(
            value[0],value[1],value[2]
        )
    )
    A2_windows=value[6].split("*")
    D5_windows=value[7].split("*")
    At_windows=value[8].split("*")
    Dt_windows=value[9].split("*")
    #* mosaic区域内的基因
    Dt_genes=Dt_synGene.loc[ 
            Dt_synGene['HC04'].isin(
                Dt_windows
            )
    ].iloc[:,36:]
    At_genes=At_synGene.loc[ 
            At_synGene['HC04'].isin(
                At_windows
            )
    ].iloc[:,36:]
    A2_genes=A2_synGene.loc[ 
        A2_synGene['J85'].isin(
            A2_windows
        )
    ].iloc[:,16:]
    D5_genes=D5_synGene.loc[
        D5_synGene['D5'].isin(
            D5_windows
        )
    ].iloc[:,2:]
    #! mosaic区域内为空?
    #* 多个windows中相同的OrthoGroup 进行合并
    D5_gene_merge=mergeOrthoGroupData(
        D5_genes,D5_genes.columns
    )
    A2_gene_merge=mergeOrthoGroupData(
        A2_genes,A2_genes.columns
    )
    At_newColumns=[i if i =="orthId" else "{}_At".format(i) for i in At_genes.columns ]
    Dt_newColumns=[i if i =="orthId" else "{}_Dt".format(i) for i in Dt_genes.columns ]
    At_gene_merge=mergeOrthoGroupData(
        At_genes,At_newColumns
    )
    Dt_gene_merge=mergeOrthoGroupData(
        Dt_genes,Dt_newColumns
    )
    mergeData1=pd.merge(
    A2_gene_merge,
    D5_gene_merge,
    on=['orthId'],
    how='outer'
    )
    mergeData2=pd.merge(
        mergeData1,
        At_gene_merge,
        on=['orthId'],
        how='outer'
    )
    #* 15+D5+35_At+35_Dt
    mergeData3=pd.merge(
        mergeData2,
        Dt_gene_merge,
        on=['orthId'],
        how='outer'
    ) 
    #* 调整merge后的columns顺序
    mergeData3=mergeData3[ 
        [ 
            'orthId', 'J85_gene', 'DC001_gene', 'DC053_gene',
            'DC086_gene', 'DC089_gene', 'DC097_gene', 'DC113_gene', 'DC119_gene',
            'DC133_gene', 'DC146_gene', 'DC151_gene', 'DC165_gene', 'DC175_gene',
            'DC212_gene', 'J98_gene', 'D5_gene', 'HC04_gene_At', 'HC15_gene_At',
            'HW03_gene_At', 'HW05_gene_At', 'HW06_gene_At', 'HW07_gene_At',
            'P01_gene_At', 'P02_gene_At', 'P04_gene_At', 'P19_gene_At',
            'P20_gene_At', 'TW007_gene_At', 'TW013_gene_At', 'TW026_gene_At',
            'TW029_gene_At', 'TW031_gene_At', 'TW055_gene_At', 'TW064_gene_At',
            'TW075_gene_At', 'TW077_gene_At', 'TW091_gene_At', 'TW094_gene_At',
            'TW100_gene_At', 'TW134_gene_At', 'XJ74_gene_At', 'XZ142_gene_At',
            'ZY006_gene_At', 'ZY10_gene_At', 'ZY184_gene_At', 'ZY236_gene_At',
            'ZY238_gene_At', 'ZY354_gene_At', 'ZY381_gene_At', 'ZY384_gene_At',
            'ZY461_gene_At', 'HC04_gene_Dt', 'HC15_gene_Dt', 'HW03_gene_Dt',
            'HW05_gene_Dt', 'HW06_gene_Dt', 'HW07_gene_Dt', 'P01_gene_Dt',
            'P02_gene_Dt', 'P04_gene_Dt', 'P19_gene_Dt', 'P20_gene_Dt',
            'TW007_gene_Dt', 'TW013_gene_Dt', 'TW026_gene_Dt', 'TW029_gene_Dt',
            'TW031_gene_Dt', 'TW055_gene_Dt', 'TW064_gene_Dt', 'TW075_gene_Dt',
            'TW077_gene_Dt', 'TW091_gene_Dt', 'TW094_gene_Dt', 'TW100_gene_Dt',
            'TW134_gene_Dt', 'XJ74_gene_Dt', 'XZ142_gene_Dt', 'ZY006_gene_Dt',
            'ZY10_gene_Dt', 'ZY184_gene_Dt', 'ZY236_gene_Dt', 'ZY238_gene_Dt',
            'ZY354_gene_Dt', 'ZY381_gene_Dt', 'ZY384_gene_Dt', 'ZY461_gene_Dt'
        ]
    ]
    for val in mergeData3.values:
        #* mosaic区域内的基因信息
        AtGenes=val[17:52]
        DtGenes=val[52:87]
        A2Genes=val[1:16]
        D5Genes=val[16:17]
        AtGeneCount=len([i for i in AtGenes if not i=="." and not pd.isna(i)])
        DtGeneCount=len([i for i in DtGenes if not i=="." and not pd.isna(i)])
        A2GeneCount=len([i for i in A2Genes if not i=="." and not pd.isna(i)])
        D5GeneCount=len([i for i in D5Genes if not i=="." and not pd.isna(i)])
        GeneData.append(
            ( 
                ["{}:{}-{}".format(value[0],value[1],value[2]),"{}:{}-{}".format(value[3],value[4],value[5]),AtGeneCount,DtGeneCount,A2GeneCount,D5GeneCount]+list(val)
            )
        )
#*Pan基因在四个亚组中多少份材料中存在
GeneData=pd.DataFrame(
        GeneData,columns=[ 
          "A_region","D_region",'AtGeneCount','DtGeneCount','A2GeneCount','D5GeneCount'
        ]+list(mergeData3.columns)
    )

GeneData.to_csv(
    sys.argv[2],header=True,index=False,sep="\t",
    na_rep="."
)