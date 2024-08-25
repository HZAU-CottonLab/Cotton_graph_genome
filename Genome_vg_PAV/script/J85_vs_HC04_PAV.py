'''
Descripttion:  获取At和A2之间的PAV, 存在跨染色体的情况
version: 
Author: zpliu
Date: 2023-06-27 11:39:15
LastEditors: zpliu
LastEditTime: 2023-07-18 15:49:03
@param: 
'''
import pandas as pd 
import re 
import sys 
bpMatrix=pd.read_csv(sys.argv[1],header=0,index_col=None,sep="\t")
searchPattern=r'([^,]+),([^:]+):(-?[^-]+)-(-?[0-9NA]+)'
#* 获取J85与HC04在哪些区域发生PAVs
#! J85的Path同时与HC04的另外两个path存在SVs?
#* 保留基因型为0的path，否则丢弃掉
#! 例如匹配染色体位于10列，则使用10
matchPathColumns=[int(i)-1 for i in sys.argv[2].split(",")]

out=[]
for value in bpMatrix.values:
    J85_SV_Length=value[4]
    J85Chr,J85Start,J85End=value[0:3]
    MatchData={}
    for matchColumn in matchPathColumns:
        #* 存在染色体易位，筛选相对应的染色体发生SVs区段
        genotype,Chr,Site1,Site2=re.findall(
                searchPattern,value[matchColumn]
            )[0]
        MatchData[Chr]=MatchData.get(Chr,[genotype,Site1,Site2])
    #! 筛选对应的染色体
    GenotypeArray=[matchItem[0] for matchItem in MatchData.values()]
    if len(list(set(GenotypeArray)))==1 and '.' in GenotypeArray:
        #* 在其他Path中都是.,可能不在所选择的Path内
        continue
    elif '0' in GenotypeArray:
        #* 存在相同的序列
        matchGenotype='0'
    elif len(list(set(GenotypeArray)))==2 and '.' in GenotypeArray:
        matchGenotype=[i for i in GenotypeArray if i!='.'][0]
    elif len(list(set(GenotypeArray)))==1 and '.' not in GenotypeArray:
        matchGenotype=GenotypeArray[0]
    else:
        #* 在多个path上都有映射,用第一个!
        matchGenotype=GenotypeArray[0]

    for Chr in MatchData.keys():
        if  MatchData[Chr][0]==matchGenotype:
            HC04Chr=Chr
            HC04_genotype,HC04Site1,HC04Site2=MatchData[Chr]
    if HC04_genotype=='0':
        HC04_SV_Length= value[4]
    else:
        HC04_SV_Length=int(value[5].split(",")[int(HC04_genotype)-1]) 
    #* HC04坐标
    HC04Start=min(int(HC04Site1),int(HC04Site2))
    HC04End=max(int(HC04Site1),int(HC04Site2))
    #! 如果SVs的序列长度差不多，则任务HC04和J85只有SNPs/Indel的差异
    if abs(J85_SV_Length-HC04_SV_Length)<50:
        pass 
    else:
        out.append(
            (
                J85Chr,J85Start,J85End,HC04Chr,HC04Start,HC04End,value[3],J85_SV_Length,HC04_SV_Length,value[6]
            )
        )
out=pd.DataFrame(out)
out.to_csv(sys.argv[3],header=False,index=False,sep="\t")
