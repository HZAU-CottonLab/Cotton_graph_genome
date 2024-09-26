###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-08-07 09:19:31
 # @LastEditors: zpliu
 # @LastEditTime: 2023-08-07 15:27:49
 # @@param: 
### 

#* 获得指定两个基因组之间的PSL文件
halLiftover='/cotton/Liuzhenping/Pan-genome/software/cactus-bin-v2.5.1/bin/halLiftover'
halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/A2_genome/panGenome-Chr01/Chr01.full.hal'
echo "test" | awk '{print "A01\t114002\t136420" }' |
    hal2maf ${halFile} --refGenome HC04-A01 --targetGenomes HC15-A01 --refTargets stdin test.maf

#* 基于halLiftover 提取HC04窗口在其他基因组中的坐标
halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/A_genome_V2/panGenome-A01/A01.full.hal'
echo "test" | awk '{print "A01\t1579500\t1602851" }' |
    $halLiftover --outPSL ${halFile} HC04-A01 stdin P04-A01 test.bed 


#* A2
halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/A2_genome/panGenome-Chr01/Chr01.full.hal'
echo "test" | awk '{print "Chr01\t949136\t954499" }' |
    $halLiftover --outPSL ${halFile} J85-Chr01 stdin DC001-Chr01 test.bed 
  