'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-10-21 17:30:37
LastEditors: zpliu
LastEditTime: 2023-10-21 19:47:38
@param: 
'''
'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-10-21 17:29:53
LastEditors: zpliu
LastEditTime: 2023-10-21 17:32:58
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

#*A2和At的所有片段
A2_At_all_windows=pd.read_csv(
    "A2_At_AllwindowAnnotate.txt",
    header=None,index_col=None,sep="\t"
)
#*A2和At的所有片段
D5_Dt_all_windows=pd.read_csv(
    "D5_Dt_AllwindowAnnotate.txt",
    header=None,index_col=None,sep="\t"
)


RegionMatch=[]
regionMatchDict={
    "":"."
}

At_Dt_SYN=pd.read_csv(
    sys.argv[1],header=None,index_col=None,sep="\t"
)
for value in At_Dt_SYN.values:
    A2_Chr,A2_start,A2_end=value[0:3]
    D5_Chr,D5_start,D5_end=value[3:6]
    #* 基于A2和D5之间的坐标筛选区间
    A2_query=pd.DataFrame(
        [(A2_Chr,A2_start,A2_end)]
    )
    D5_query=pd.DataFrame(
        [(D5_Chr,D5_start,D5_end)]
    )
    #! 基于At和Dt之间的坐标进行筛选对于的SYN和hyper区域
    A2_data=intersectBed(
        A2_query,
        A2_At_all_windows.loc[A2_At_all_windows[3]!="."][ 
                [3,4,5,0,1,2,6]
            ]
    )
    D5_data=intersectBed(
        D5_query,
        D5_Dt_all_windows.loc[ 
            D5_Dt_all_windows[3]!="."
        ][ 
            [3,4,5,0,1,2,6]
        ]
    )
    #* 获取对应的，比较模糊的共线性区块
    Dt_region="*".join(
        list(
            D5_data.loc[ 
            D5_data[3]!="."
        ][[3,4,5]].apply(
            lambda x:"{}:{}-{}".format(
                x[3],x[4],x[5]
            ),axis=1
        ).values
        )
    )
    D5_region="*".join(
        list(
            D5_data.loc[ 
            D5_data[6]!="."
        ][[6,7,8]].apply(
            lambda x:"{}:{}-{}".format(
                x[6],x[7],x[8]
            ),axis=1
        ).values
        )
    )
    At_region="*".join(
        list(
            A2_data.loc[ 
            A2_data[3]!="."
        ][[3,4,5]].apply(
            lambda x:"{}:{}-{}".format(
                x[3],x[4],x[5]
            ),axis=1
        ).values
        )
    )
    A2_region="*".join(
        list(
            A2_data.loc[ 
            A2_data[6]!="."
        ][[6,7,8]].apply(
            lambda x:"{}:{}-{}".format(
                x[6],x[7],x[8]
            ),axis=1
        ).values
        )
    )
    At_region=regionMatchDict.get(
        At_region,At_region
    )
    A2_region=regionMatchDict.get(
        A2_region,A2_region
    )
    Dt_region=regionMatchDict.get(
        Dt_region,Dt_region
    )
    D5_region=regionMatchDict.get(
        D5_region,D5_region
    )
    #! 对相应区间内的所有基因进行overlap
    RegionMatch.append( 
        ( A2_Chr,A2_start,A2_end,D5_Chr,D5_start,D5_end,A2_region,D5_region,At_region,Dt_region)
    )

RegionMatch=pd.DataFrame(
    RegionMatch,
    columns=[ 
        'A2Chr','A2Start','A2End','D5Chr','D5Start','D5End',
        'A2_region','D5_region','At_region','Dt_region'
    ]
)
RegionMatch.to_csv(
    sys.argv[2],header=True,index=False,sep="\t"
)