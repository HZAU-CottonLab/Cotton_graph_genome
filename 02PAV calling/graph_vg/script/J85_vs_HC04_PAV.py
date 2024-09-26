'''
Descripttion:  Get the PAV between At and A2, which is cross-chromosome.
version: 
Author: zpliu
Date: 2023-06-27 11:39:15
LastEditors: zpliu
LastEditTime: 2024-09-26 09:46:17
@param: 
'''
import pandas as pd 
import re 
import sys 
bpMatrix=pd.read_csv(sys.argv[1],header=0,index_col=None,sep="\t")
searchPattern=r'([^,]+),([^:]+):(-?[^-]+)-(-?[0-9NA]+)'
#* Obtain the areas where PAVs of J85 and HC04 occur.
#* 保留基因型为0的path，否则丢弃掉
matchPathColumns=[int(i)-1 for i in sys.argv[2].split(",")]

out=[]
for value in bpMatrix.values:
    J85_SV_Length=value[4]
    J85Chr,J85Start,J85End=value[0:3]
    MatchData={}
    for matchColumn in matchPathColumns:
        # * Chromosome translocation exists, and the corresponding chromosome SVs segment is screened.
        genotype,Chr,Site1,Site2=re.findall(
                searchPattern,value[matchColumn]
            )[0]
        MatchData[Chr]=MatchData.get(Chr,[genotype,Site1,Site2])
    GenotypeArray=[matchItem[0] for matchItem in MatchData.values()]
    if len(list(set(GenotypeArray)))==1 and '.' in GenotypeArray:
        continue
    elif '0' in GenotypeArray:
        matchGenotype='0'
    elif len(list(set(GenotypeArray)))==2 and '.' in GenotypeArray:
        matchGenotype=[i for i in GenotypeArray if i!='.'][0]
    elif len(list(set(GenotypeArray)))==1 and '.' not in GenotypeArray:
        matchGenotype=GenotypeArray[0]
    else:
    
        matchGenotype=GenotypeArray[0]

    for Chr in MatchData.keys():
        if  MatchData[Chr][0]==matchGenotype:
            HC04Chr=Chr
            HC04_genotype,HC04Site1,HC04Site2=MatchData[Chr]
    if HC04_genotype=='0':
        HC04_SV_Length= value[4]
    else:
        HC04_SV_Length=int(value[5].split(",")[int(HC04_genotype)-1]) 
    #* HC04
    HC04Start=min(int(HC04Site1),int(HC04Site2))
    HC04End=max(int(HC04Site1),int(HC04Site2))
   
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
