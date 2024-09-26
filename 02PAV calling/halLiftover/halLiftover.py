'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-07-26 19:45:12
LastEditors: zpliu
LastEditTime: 2023-08-04 14:52:06
@param: 
'''
from tempfile import NamedTemporaryFile
import pandas as pd 
import subprocess
import os 

def halLiftover(queryBed,halFile,queryGenome,targetGenome):
    '''
    @queryBed: pd.DataFrame
    @halFile:
    @queryGenome:
    @targetGenome: 
    '''
    halbin='/cotton/Liuzhenping/Pan-genome/software/cactus-bin-v2.5.1/bin/halLiftover'
    queryBedFile=NamedTemporaryFile(mode="w+")
    targetFile=NamedTemporaryFile(mode="w+")
    queryBed.to_csv(
        queryBedFile.name,header=False,index=False,sep="\t"
    )
    #* 使用halLiftover获取坐标
    commond='{} --outPSL {} {} {} {} {} '.format(
        halbin,halFile,queryGenome,queryBedFile.name,targetGenome,targetFile.name
    )
    subprocess.run(
       args=[commond],
        check=True,
        shell=True
    )
    #* 获取其所在的区间
    if os.path.getsize(targetFile.name)>0:
        PSLout=pd.read_csv(
            targetFile.name,
            header=None,index_col=None,sep="\t"
        )
    else:
        raise Exception("query Chrom not in hal or can't liftover this region")
    queryBedFile.close()
    targetFile.close()
    return PSLout

if __name__ =="__main__":
    #* 查询指定区域在另外一个基因组的liftover坐标
    queryBed=pd.DataFrame(
        [ 
            ("A13",99781907,99781914)
        ]
    )
    halLiftover(
    queryBed,'/public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/A_genome_V2/panGenome-A13/A13.full.hal',
    "HC04-A13","TW100-A13"
    )