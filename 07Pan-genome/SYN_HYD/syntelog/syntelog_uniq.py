'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-10-23 17:30:58
LastEditors: zpliu
LastEditTime: 2023-10-23 17:40:53
@param: 
'''
#TODO 基于三个syntelogs数据集，进行去重
from operator import index
import pandas as pd 
import sys 
import networkx as nx
from itertools import combinations
import re 

#! 按照OrthoID进行处理
def getVariantCluster(NodeData):
    MapDic={}
    for item in NodeData:
        SNPid1,SNPid2=item
        MapDic[SNPid1]=MapDic.get(SNPid1,[SNPid2])
        MapDic[SNPid2]=MapDic.get(SNPid1,[SNPid1])
    #* 获得所有的Node list
    SNPidlist=MapDic.keys()
    #TODO 基于连通图，提取子节点
    SNPcluster=[]
    G = nx.Graph()
    for node in SNPidlist:
        G.add_node(node)
    for link in NodeData:
        G.add_edge(link[0], link[1])
    for c in nx.connected_components(G):
        #* 获得每个subgraph内的节点信息
        nodeSet=G.subgraph(c).nodes()
        #* 存在相互联系的节点列表
        SNPcluster.append(nodeSet)
    return SNPcluster

def Order_sample_List(geneNodes,OrthId):
    '''
    #* 按照顺序排列
    '''
    SampleOrder=[ 
            'J85', 'DC001', 'DC053',
            'DC086', 'DC089', 'DC097', 'DC113', 'DC119',
            'DC133', 'DC146', 'DC151', 'DC165', 'DC175',
            'DC212', 'J98', 'Grai', 'HC04_A', 'HC15_A',
            'HW03_A', 'HW05_A', 'HW06_A', 'HW07_A',
            'P01_A', 'P02_A', 'P04_A', 'P19_A',
            'P20_A', 'TW007_A', 'TW013_A', 'TW026_A',
            'TW029_A', 'TW031_A', 'TW055_A', 'TW064_A',
            'TW075_A', 'TW077_A', 'TW091_A', 'TW094_A',
            'TW100_A', 'TW134_A', 'XJ74_A', 'XZ142_A',
            'ZY006_A', 'ZY10_A', 'ZY184_A', 'ZY236_A',
            'ZY238_A', 'ZY354_A', 'ZY381_A', 'ZY384_A',
            'ZY461_A', 'HC04_D', 'HC15_D', 'HW03_D',
            'HW05_D', 'HW06_D', 'HW07_D', 'P01_D',
            'P02_D', 'P04_D', 'P19_D', 'P20_D',
            'TW007_D', 'TW013_D', 'TW026_D', 'TW029_D',
            'TW031_D', 'TW055_D', 'TW064_D', 'TW075_D',
            'TW077_D', 'TW091_D', 'TW094_D', 'TW100_D',
            'TW134_D', 'XJ74_D', 'XZ142_D', 'ZY006_D',
            'ZY10_D', 'ZY184_D', 'ZY236_D', 'ZY238_D',
            'ZY354_D', 'ZY381_D', 'ZY384_D', 'ZY461_D'
        ]
    reorderGeneList=[]
    for sample in SampleOrder:
        searchPattern=r'{}.*'.format(sample)
        matchData=[re.findall(searchPattern,genes)[0] for genes in geneNodes if re.findall(searchPattern,genes)]
        if len(matchData)==0:
            reorderGeneList.append(".")
        else:
            reorderGeneList.append(
                "-".join(matchData)
            )
    A2_count=len([i for i in reorderGeneList[0:15] if i!="."])
    D5_count=len([i for i in reorderGeneList[15:16] if i!="."])
    At_count=len([i for i in reorderGeneList[16:51] if i!="."])
    Dt_count=len([i for i in reorderGeneList[51:86] if i!="."])
    return [OrthId,At_count,Dt_count,A2_count,D5_count]+reorderGeneList

mergeData=pd.read_csv(
    "./network_inputData.txt",
    header=0,index_col=None,sep="\t"
)
queryOrthoList=pd.read_csv(
    sys.argv[1],header=None,index_col=None,sep="\t"
)

All_syntelog_mergeData=[]
for queryOrthoId in queryOrthoList[0].values:
    #* 构造所有的links
    syntelogList=[]
    for val in mergeData.loc[ 
                        mergeData['orthId']==queryOrthoId
                    ].values:
        sample_with_geneList=[i for i in val[5:] if i!="."]
        if len(sample_with_geneList)==1:
            #! 该共线性区域，只在一份材料中注释出了Gene，但是该基因在整个基因组水平存在同源基因
            pass 
        for A,B in combinations(sample_with_geneList,2):
            syntelogList.append(
                (A,B)
            )
    #* 合并后的pan-syntelogs
    unique_syntelogList=getVariantCluster(syntelogList)
    #! 按照样本顺序进行排列后的结果
    for syntelogs in unique_syntelogList:
        #! 不同材料的syntelog node进行排序
        orderSyntelog=Order_sample_List( 
            syntelogs,queryOrthoId
        )
        All_syntelog_mergeData.append(
            orderSyntelog
        )

All_syntelog_mergeData=pd.DataFrame(
    All_syntelog_mergeData,
    columns= [ 
            'orthId','AtGeneCount','DtGeneCount','A2GeneCount','D5GeneCount',
             'J85_gene', 'DC001_gene', 'DC053_gene',
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
        ])

All_syntelog_mergeData.to_csv(sys.argv[2],header=True,index=False,sep="\t")
