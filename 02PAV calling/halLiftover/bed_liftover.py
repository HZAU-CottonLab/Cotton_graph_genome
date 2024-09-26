'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-07-27 10:46:07
LastEditors: zpliu
LastEditTime: 2023-09-14 16:34:31
@param: 
'''

import sys
sys.path.append("/public/home/zpliu/Pan-genome/SV_parallele_V2/GenomeCompare")
sys.path.append("/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/")
import pandas as pd
from divergenceRegion.bedtools import intersectBed
from halLiftover import halLiftover
import logging
logging.basicConfig(
    level=logging.INFO
)
def bed_liftover(halFile:str,queryBed:pd.DataFrame,queryGenome:str,targetGenome:str):
    try:
        logging.info(
            "start to liftover..."
        )
        liftover_PSL = halLiftover(
            queryBed, halFile, queryGenome, targetGenome
        )
    except Exception as e:
        logging.info(e)
        return False
    intersectedData = intersectBed(
        queryBed, liftover_PSL[[9, 11, 12, 13, 15, 16, 0, 2, 8]]
    )
    # * 过滤query Item存在交集的坐标
    MatchData = []
    for window, windowLiftover in intersectedData.groupby([0, 1, 2]):
        if windowLiftover.iloc[0, 3] == ".":
            MatchData.append(
                (
                    window[0], int(window[1]), int(window[2]),
                    window[0], int(window[1]), int(window[2]), ".", -1, -1, 0, "."
                )
            )
        else:
            windowLiftover[12] = windowLiftover.apply(
                lambda x: int(x[9])+int(x[10]), axis=1
            )
            windowLiftover = windowLiftover.sort_values(by=[12], ascending=False)
            # * 过滤掉query小于intersect的query, 防止不同的query windows存在交集
            windowLiftover = windowLiftover.loc[
                windowLiftover.apply(
                    lambda x:True if int(x[1]) <= int(x[4]) and int(x[2]) >= int(x[5]) else False, axis=1
                )
            ]
            if windowLiftover.empty:
                MatchData.append(
                    (
                        window[0], int(window[1]), int(window[2]), window[0], int(
                            window[1]), int(window[2]), ".", -1, -1, 0, "."
                    )
                )
                continue
            # * 每个windows在该基因组中映射的坐标
            for val in windowLiftover.values:
                MatchData.append(
                    (
                        window[0], int(window[1]), int(window[2]),
                        val[3], val[4], val[5],val[6],val[7],val[8],
                        val[12],val[11]
                    )
                )
    return pd.DataFrame(MatchData, columns=['queryChr', "queryStart", "queryEnd",
                                                'MatchChr', 'MatchStart', 'MatchEnd',
                                                'targetChr', 'targetStart', 'targetEnd', 'MatchCount', 'Stand'])
if __name__=="__main__" :
    # * 文件包含3列
    halFile = sys.argv[1]
    queryBed = pd.read_csv(
        sys.argv[2], header=None, index_col=None, sep="\t"
    )
    queryGenome = sys.argv[3]
    targetGenome = sys.argv[4]
    outFile = sys.argv[5]

    MatchData = bed_liftover(
        halFile,queryBed,queryGenome,targetGenome
    )
    MatchData.to_csv(
        outFile, header=True, index=False, sep="\t"
    )

    logging.info(
        "Completed!"
    )
    # * 将query的Bed文件与liftover文件取交集
