import pandas as pd
core_num = pd.read_table("core_num.txt",header=None)
core_num["mean"] = core_num.iloc[:,1:].mean(axis = 1)
core_num["STD"] =  core_num.iloc[:,1:].std(axis = 1)
core_num = core_num.loc[:,[0,"mean","STD"]]

pan_num = pd.read_table("pan_num.txt",header= None)
pan_num["mean"] = pan_num.iloc[:,1:].mean(axis = 1)
pan_num["STD"] = pan_num.iloc[:,1:].std(axis = 1)
pan_num = pan_num.loc[:,[0,"mean","STD"]]

core_num["type"] = "core"
pan_num["type"] = "pan"
merge = pd.concat([core_num,pan_num])
merge.to_csv("core_pan.txt",index=False,sep="\t")
