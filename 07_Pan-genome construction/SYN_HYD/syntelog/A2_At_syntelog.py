'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-10-23 11:58:13
LastEditors: zpliu
LastEditTime: 2023-10-23 15:08:35
@param: 
'''
import pandas as pd 
import sys 
import logging
logging.basicConfig(
    level=logging.INFO
)

A2_At_allWindows=pd.read_csv(
    sys.argv[1],
    header=None,index_col=None,sep="\t"
)
A2_syntelogs=pd.read_csv(
    "./A2_all_window_liftover_15Gene.txt",
    header=0,index_col=None,sep="\t"
)
At_syntelogs=pd.read_csv(
    "./At_all_window_liftover_35Gene.txt",
    header=0,index_col=None,sep="\t"
)



A2_At_syntelogs_list=[]
for val in A2_At_allWindows.values:
    A2_region="{}:{}-{}".format( 
        val[0],val[1],val[2]
    )
    At_region="{}:{}-{}".format(
        val[3],val[4],val[5]
    )
    #* 
    logging.info(
        "{}-{}>>>>>".format(
            A2_region,At_region
        )
    )
    A2_query_syntelogs=A2_syntelogs.loc[ 
                A2_syntelogs['J85']==A2_region
            ].iloc[:,16:].drop_duplicates()
    At_query_syntelogs=At_syntelogs.loc[
            At_syntelogs['HC04']==At_region
        ].iloc[:,36:].drop_duplicates()
    if A2_query_syntelogs.empty and At_query_syntelogs.empty:
        continue 
    elif A2_query_syntelogs.empty: 
        for val in At_query_syntelogs.values:
            AtgeneCount=len([i for i in val[1:36] if i!="." and not pd.isna(i)])
            A2_At_syntelogs_list.append(
                [A2_region,At_region,AtgeneCount,0,0,0]+[val[0]]+["."]*16+list(val[1:36])+["."]*35
            )
    elif At_query_syntelogs.empty:
        #! hyper区域默认没有syntelogs
        for val in A2_query_syntelogs.values:
            A2geneCount=len([i for i in val[1:16] if i!="." and not pd.isna(i)])
            A2_At_syntelogs_list.append(
                [A2_region,At_region,0,0,A2geneCount,0]+list(val[0:16])+["."]+["."]*70
            ) 
    else:
        merge_syntelogs=pd.merge(
            A2_query_syntelogs,
            At_query_syntelogs,
            left_on=['orthId'],
            right_on=['orthId'],
            how='outer'
        )
        merge_syntelogs=merge_syntelogs[ 
            [ 
            'orthId', 'J85_gene', 'DC001_gene', 'DC053_gene',
            'DC086_gene', 'DC089_gene', 'DC097_gene', 'DC113_gene', 'DC119_gene',
            'DC133_gene', 'DC146_gene', 'DC151_gene', 'DC165_gene', 'DC175_gene',
            'DC212_gene', 'J98_gene', 'HC04_gene', 'HC15_gene',
            'HW03_gene', 'HW05_gene', 'HW06_gene', 'HW07_gene',
            'P01_gene', 'P02_gene', 'P04_gene', 'P19_gene',
            'P20_gene', 'TW007_gene', 'TW013_gene', 'TW026_gene',
            'TW029_gene', 'TW031_gene', 'TW055_gene', 'TW064_gene',
            'TW075_gene', 'TW077_gene', 'TW091_gene', 'TW094_gene',
            'TW100_gene', 'TW134_gene', 'XJ74_gene', 'XZ142_gene',
            'ZY006_gene', 'ZY10_gene', 'ZY184_gene', 'ZY236_gene',
            'ZY238_gene', 'ZY354_gene', 'ZY381_gene', 'ZY384_gene',
            'ZY461_gene'
            ]
        ]
        for val in merge_syntelogs.values:
            AtgeneCount=len([i for i in val[16:] if i!="." and not pd.isna(i)])
            A2geneCount=len([i for i in val[1:16] if i!="." and not pd.isna(i)])
            A2_At_syntelogs_list.append(
                [A2_region,At_region,AtgeneCount,0,A2geneCount,0]+list(val[0:16])+["."]+list(val[16:])+["."]*35
            )
A2_At_syntelogs_list=pd.DataFrame(
    A2_At_syntelogs_list,
    columns=[ 
           "A_region","D_region",'AtGeneCount','DtGeneCount','A2GeneCount','D5GeneCount',
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
    )
A2_At_syntelogs_list.to_csv(
    sys.argv[2],header=True,index=False,sep="\t",na_rep="."
)