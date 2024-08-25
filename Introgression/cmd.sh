###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-11-05 10:30:45
 # @LastEditors: zpliu
 # @LastEditTime: 2023-11-22 22:37:24
 # @@param: 
### 
#! lijianyin Genome Biolog QTL  TM1.0坐标转HC04坐标
#? 
script="/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/bed_liftover.py"
halFile="/public/home/zpliu/Pan-genome/Cactus-Pan/Allopolyploids/genomeAlignment-pg/genomeAlignment-pg.full.hal"
queryGenome="TM1_V1"
targetGenome="HC04"
queryBedFile="./TM1_V1.0_jyli_GB.txt"
outFile="HC04_jyli_GB.txt"
python ${script} ${halFile} ${queryBedFile} ${queryGenome} ${targetGenome} ${outFile}

#* 用黄越凡附表中的QTL进行了相应的分析
script="/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/bed_liftover.py"
halFile="/public/home/zpliu/Pan-genome/Cactus-Pan/Allopolyploids/genomeAlignment-pg/genomeAlignment-pg.full.hal"
queryGenome="TM1_V2"
targetGenome="HC04"
queryBedFile="./11"
outFile="HC04_YFhuang.txt"
python ${script} ${halFile} ${queryBedFile} ${queryGenome} ${targetGenome} ${outFile}





#* 渐渗结果所在的目录：
#? /cotton/JianyingLi/PAN_Graph/AD1/06.Dom_Imp/PanSV_IS_V1 


#* 三种渐渗区间中，所包含的hyper序列的占比是否有差异(没有差异)



#* 将渐渗区间与表型GWAS的QTL取交集
#? 黄越凡 搜集QTLs
/public/home/jqyou/out/pan_genome_out/supplementary_data/HC04_published_QTL.txt
#? NG 2023 搜集QTLs
/public/home/jqyou/out/pan_genome_out/supplementary_data/hg_QTL_new_HC04_V1.txt



#* 对渐渗区域内是位于geneic还是非geneic区域进行注释
#* 渐渗区间与其他区间相比，基因的PAV频率偏差即在野生棉和栽培棉中基因PAV的频率不一样？
#! gene PAV biases 
#? 与渐渗片段存在交集的gene,其在野生棉和栽培棉中的PAV biases差异越大 



#TODO 亚基因组在渐渗区间对应的比较
#? 如何获取At和Dt 1对1的保守渐渗区间；



#TODO Dt在整个基因组水平的渐渗密度高于At？
