'''
Descripttion: 基于PAV文件,计算不同材料在指定区域的多态性
version: 
Author: zpliu
Date: 2023-08-23 14:55:05
LastEditors: zpliu
LastEditTime: 2023-08-28 09:22:22
@param: 
'''
import pandas as pd 
import sys
import logging
logging.basicConfig(
    level=logging.INFO
)
#* 导入所需要的包
#* 导入所需要的包
sys.path.append("/public/home/zpliu/Pan-genome/SV_parallele_V2/GenomeCompare")
sys.path.append("/public/home/zpliu/Pan-genome/Cactus-Pan/")
from divergenceRegion.bedtools import intersectBed
from halLiftover.halLiftover import halLiftover
from itertools import combinations


SV_info=pd.read_csv(
        "/public/home/zpliu/Pan-genome/SV_parallele_V2/superPan/A2_genotype_V2/A2_SV_merge_info.txt",
    header=0,index_col=None,sep="\t"
    )
AD1_genome=pd.read_csv(
    "/public/home/zpliu/Pan-genome/SV_parallele_V2/superPan/A2_longRead_call/15.Assembly.txt",
    header=None,index_col=None,sep="\t"
)

#* 根据材料的基因型信息填充序列长度信息
imputedSV=[]
for val in SV_info.values:
    seqLength=[]
    seqMap=[val[4]]
    for i in val[5].split(","):
        #* alt seq length
        seqMap.append(
            int(i)
        )
    for genotype in val[7:]:
        if genotype=="0"  or genotype==0 or genotype=='.':
            seqLength.append(
                seqMap[0]
            ) 
        else:
            seqLength.append(
                seqMap[int(genotype)]
            )
    imputedSV.append(
        list(val[0:7])+seqLength
    )
imputedSV=pd.DataFrame(imputedSV,columns=SV_info.columns) 

#* 测试前10个windows
queryBedFile=pd.read_csv(
   sys.argv[1],header=None,index_col=None,sep="\t"
)



intersectedData=intersectBed(
    queryBedFile,imputedSV
)
intersectedData.columns=[
    'windowChr','windowStart','windowEnd',
    '#Chr', 'start', 'end', 'ID', 'refLen', 'AltLen', 'SV_level', 'J85',
    'DC001', 'DC053', 'DC086', 'DC089', 'DC097', 'DC113', 'DC119', 'DC133', 'DC146',
    'DC151', 'DC165', 'DC175', 'DC212', 'J98'
]

intersectedData=intersectedData.astype(
    {
        "windowStart":int,
        "windowEnd":int
    }
)
#! 计算windows在多个材料间的多态性

Pi_result=[]
for window,windowData in intersectedData.groupby(['windowChr','windowStart','windowEnd']):
    if windowData.iloc[0,3]==".":
        #* 该区域没有SVs的信息，默认该区域的PAV多态性为0
        Pi_result.append(
            (window[0],window[1],window[2],0)
        )
        continue
    #* 计算HC04 windows在其他基因组中对应的windows
    #! 如果该windows与SV坐标没有交集，对应的多态性值为0
    logging.info("parseing {}".format(window))
    queryChr=window[0][4:]
    halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/A2_genome/panGenome-{}/{}.full.hal'.format(
        queryChr,queryChr
    )
    #! 获得HC04 windows在其他基因组中映射的坐标
    genomeLength=[ 
      window[2]-window[1]  
    ]
    for genome in AD1_genome.iloc[1:,0]:
        #* 定位该windows在其他基因组中的坐标
        queryGenome="J85-{}".format(queryChr)
        targetGenome="{}-{}".format(genome,queryChr)
        genomeData=[]
        try:
            #* liftover 该windows在其他基因组中的坐标
            liftoverData=halLiftover(
                pd.DataFrame([(queryChr,window[1],window[2])]),halFile,queryGenome,targetGenome
            )
        except Exception as e:
            
            genomeLength.append(
                0
            )
            continue
        if liftoverData.empty:
            genomeLength.append(
                    0
                )
        else:
            #* liftover 的windows大小
            genomeLength.append(
                    int(liftoverData.iloc[0,16])-int(liftoverData.iloc[0,15])
                )
    sampleWindowLength=dict(zip(AD1_genome[0],genomeLength))
    #* 根据任意两个基因组之间PAV序列的差异占整个windows的长度，计算Pi值
    Pi=0
    for sample1,sample2 in combinations(AD1_genome[0],2):
        #* 两个材料累积的序列差异
        seqDifference=windowData[[sample1,sample2]].apply(
            lambda x:abs(int(x[0])-int(x[1])),axis=1
        ).sum()
        sample1Length=sampleWindowLength.get(sample1)
        sample2Length=sampleWindowLength.get(sample2)
        if sample1Length==0 or sample2Length==0:
            #* 完全差异的windows
            Pi+=1/(15*15)
        else:
            #* 累积的差异程度
            differentRatio1=seqDifference/sample1Length
            differentRatio2=seqDifference/sample2Length
            if differentRatio1>1:
                '''#* 对齐比较差的区域
                '''
                differentRatio1=1
            if differentRatio2>1:
                differentRatio2=1
            Pi+=(differentRatio2+differentRatio1)/(2*15*15)
    Pi_result.append(
            (window[0],window[1],window[2],Pi*2*15/14)
        )
Pi_result=pd.DataFrame(Pi_result,columns=['Chr','start','end','Pi'])
Pi_result.to_csv(
    sys.argv[2],
    header=True,index=False,sep="\t"
)