'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-07-13 15:01:36
LastEditors: zpliu
LastEditTime: 2024-09-26 10:44:37
@param: 
'''
import pandas as pd 
import pysam
import re 
from bedtools import subtractBed_genome,merge_bed,intersectBed


if __name__=="__main__":
    #根据共线性坐标，将共线性区间进行合并
    syntenyRegion=pd.read_csv(
        "J85_vs_HC04_At_synteny.txt",
        header=None,index_col=None,sep="\t"
    )
    #* 合并存交集的共线性区间
    query_synteny=syntenyRegion.sort_values(by=[0,1,2])[[0,1,2]]
    J85_mergeBed=merge_bed(query_synteny)
    #* 合并HC04的坐标
    query_synteny=syntenyRegion.sort_values(by=[3,4,5])[[3,4,5]]
    HC04_mergeBed=merge_bed(query_synteny)

    #* J85合并的共线性与未合并前的共线性坐标取交集
    J85_intersectedData=intersectBed(J85_mergeBed,syntenyRegion)
    HC04_intersectedData=intersectBed(HC04_mergeBed,syntenyRegion[[3,4,5,0,1,2,6]])

    #* 合并两个共线性的结果
    J85_intersectedData.columns=[
        "Chr_J85_merge","Start_J85_merge","End_J85_merge",
        "Chr_J85_syn","Start_J85_syn","End_J85_syn",
        "Chr_HC04_syn","Start_HC04_syn","End_HC04_syn","Stand"
    ]
    HC04_intersectedData.columns=[
        "Chr_HC04_merge","Start_HC04_merge","End_HC04_merge",
        "Chr_HC04_syn","Start_HC04_syn","End_HC04_syn",
        "Chr_J85_syn","Start_J85_syn","End_J85_syn",
        'Stand'
    ]
    #* 转变数据类型
    J85_intersectedData=J85_intersectedData.astype(
        {
            "Start_HC04_syn":int,
            "End_HC04_syn":int
        }
    )
    HC04_intersectedData=HC04_intersectedData.astype(
        {
        "Start_J85_syn" :int,
        "End_J85_syn" :int,
        }
    )
    #* 两个基因的共线性区间
    mergeData=pd.merge(
    J85_intersectedData,
    HC04_intersectedData,
        on=[
            "Chr_J85_syn","Start_J85_syn","End_J85_syn",
            "Chr_HC04_syn","Start_HC04_syn","End_HC04_syn","Stand"
        ]
    )
    #* 两个基因组各自特异性的divergence区间
    AD1_substractBed=subtractBed_genome(
    "HC04-softMasked.fa.fai",
        mergeData[['Chr_HC04_merge',"Start_HC04_merge","End_HC04_merge"]]
    )
    At_substractBed=AD1_substractBed.loc[AD1_substractBed.apply(
        lambda x: True if re.match("^HC04_A",x[0]) else False,axis=1
    )]
    J85_substractBed=subtractBed_genome(
        "J85-softMasked.fa.fai",
        mergeData[['Chr_J85_merge',"Start_J85_merge","End_J85_merge"]]
    )
    mergeData['Type']='synteny'
    J85_substractBed[3]='divergence'
    At_substractBed[3]='divergence'
    #todo 基于VG 解构的PAV关系，将两个基因组之间各自特异性的divergence区间进行合并
    PAV_list=pd.read_csv(
        "J85_vs_HC04_synteny.txt",
        header=None,sep="\t",index_col=None
    )
    #* 将对应的PAV区间与HC04基因组中所有窗口取交集
    J85_windows_intersectedData=intersectBed(
    PAV_list,J85_substractBed
    )
    J85_windows_intersectedData=J85_windows_intersectedData[
        [3,4,5,0,1,2,7,6,8,9,10,11]
    ]
    HC04_J85_window_intersectedData=intersectBed(
        J85_windows_intersectedData,At_substractBed
    )
    #* 所有的交集情况
    HC04_J85_window_intersectedData.columns=[
        'HC04_PAV_Chr','HC04_PAV_Start',"HC04_PAV_End",
        "J85_PAV_Chr","J85_PAV_Start","J85_PAV_End",
        'HC04_PAV_len',"J85_PAV_len",
        "J85_window_Chr","J85_window_Start","J85_window_End",
        "J85_windowType","HC04_window_Chr","HC04_window_Start","HC04_window_End",
        "HC04_window_Type"

    ]
    HC04_J85_window_intersectedData=HC04_J85_window_intersectedData.astype(
        {
        "HC04_PAV_len":int,
        "J85_PAV_len":int,
        "J85_window_Start":int,
        "J85_window_End":int,
        "HC04_window_Start":int,
        "HC04_window_End":int,
        }
    )
    divergenceWindowsMatch=HC04_J85_window_intersectedData.loc[ 
    (HC04_J85_window_intersectedData['J85_windowType']=="divergence")&(
        HC04_J85_window_intersectedData['HC04_window_Type']=="divergence"
    )][['J85_window_Chr',"J85_window_Start","J85_window_End",
    "HC04_window_Chr","HC04_window_Start","HC04_window_End"]].drop_duplicates()

    J85_matchWindow={}
    At_matchWindow={}
    for value in divergenceWindowsMatch.values:
        A2_divergenceId="{}:{}-{}".format(value[0],value[1],value[2])
        At_divergenceId="{}:{}-{}".format(value[3],value[4],value[5])
        #* 对应的At divergence windows
        J85_matchWindow[A2_divergenceId]=J85_matchWindow.get(
            A2_divergenceId,[]
        )
        J85_matchWindow[A2_divergenceId].append(
            (value[3],value[4],value[5])
        )
        #* 对应的At divergence windows
        At_matchWindow[At_divergenceId]=At_matchWindow.get(
            At_divergenceId,[]
        )
        At_matchWindow[At_divergenceId].append((
            value[0],value[1],value[2]
        ))
    J85_divergenceData=[]
    for divergenceSite in J85_substractBed.values:
        A2Id="{}:{}-{}".format(divergenceSite[0],divergenceSite[1],divergenceSite[2])
        mathcData=J85_matchWindow.get(A2Id)
        if mathcData and len(mathcData)==1:
            slefMap=At_matchWindow.get("{}:{}-{}".format(mathcData[0][0],mathcData[0][1],mathcData[0][2]))
            if slefMap and len(slefMap)==1 and A2Id=="{}:{}-{}".format(slefMap[0][0],slefMap[0][1],slefMap[0][2]):
                #* 互相映射
                J85_divergenceData.append(
                    [divergenceSite[0],divergenceSite[1],divergenceSite[2]]+list(mathcData[0])+["divergence"]
                )
            else:
                #*  
                J85_divergenceData.append(
                    [divergenceSite[0],divergenceSite[1],divergenceSite[2],".",-1,-1,"divergence"]
                ) 
        else:
            J85_divergenceData.append(
                [divergenceSite[0],divergenceSite[1],divergenceSite[2],".",-1,-1,"divergence"]
            ) 
    J85_divergenceData=pd.DataFrame(J85_divergenceData)

    At_divergenceData=[]
    for divergenceSite in At_substractBed.values:
        AtId="{}:{}-{}".format(divergenceSite[0],divergenceSite[1],divergenceSite[2])
        mathcData=At_matchWindow.get(AtId)
        if mathcData and len(mathcData)==1:
            slefMap=J85_matchWindow.get("{}:{}-{}".format(mathcData[0][0],mathcData[0][1],mathcData[0][2]))
            if slefMap and len(slefMap)==1 and AtId=="{}:{}-{}".format(slefMap[0][0],slefMap[0][1],slefMap[0][2]):
                At_divergenceData.append(
                    list(mathcData[0])+ [divergenceSite[0],divergenceSite[1],divergenceSite[2],"divergence"]
                )
            else:
                At_divergenceData.append(
                    [".",-1,-1,divergenceSite[0],divergenceSite[1],divergenceSite[2],"divergence"]
                ) 
        else:
            At_divergenceData.append(
                [".",-1,-1,divergenceSite[0],divergenceSite[1],divergenceSite[2],"divergence"]
            ) 
    At_divergenceData=pd.DataFrame(At_divergenceData)


    mergedivergenceData=pd.concat(
        [At_divergenceData,J85_divergenceData],axis=0
    )
    mergedivergenceData=mergedivergenceData.drop_duplicates()
    mergedivergenceData.columns=[
        'Chr_J85_merge',"Start_J85_merge","End_J85_merge",
        "Chr_HC04_merge","Start_HC04_merge","End_HC04_merge","Type"]
    #* 共线性的窗口
    syntenyMerge=mergeData[[
        'Chr_J85_merge',"Start_J85_merge","End_J85_merge",
        "Chr_HC04_merge","Start_HC04_merge","End_HC04_merge","Type"]]
    #! 最终的windows
    Finally_windows=pd.concat(
        [mergedivergenceData,syntenyMerge],axis=0
    )
    Finally_windows=Finally_windows.drop_duplicates()
    Finally_windows.to_csv(
        "J85_vs_At_all_windows.txt",
        header=True,index=False,sep="\t"
    )

    