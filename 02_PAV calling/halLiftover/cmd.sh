###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-07-27 10:54:40
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 10:15:51
# @@param:
###

#TODO Based on the bed of HC04, query its coordinates in other genome.

script='/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/bed_liftover.py'
halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/Allopolyploids/genomeAlignment-pg/genomeAlignment-pg.full.hal'
queryGenome='HC04'
targetGenome='TM1_V2'
queryBedFile='/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/test/HC04.bed'
outFile='/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/test/test.txt'

#* run
#? python with version ＞3.6
python ${script} ${halFile} ${queryBedFile} ${queryGenome} ${targetGenome} ${outFile}




script='/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/bed_liftover.py'

halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/A_genome_V2/panGenome-A01/A01.full.hal'

halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/D_genome_V2/panGenome-D01/D01.full.hal'
halFile='/public/home/zpliu/Pan-genome/Cactus-Pan/A2_genome/panGenome-Chr01/Chr01.full.hal'
queryGenome='HC04-A01'
targetGenome='P01-A01'

queryBedFile='/public/home/zpliu/Pan-genome/Cactus-Pan/halLiftover/test/HC04.bed'
outFile='test.txt'
python ${script} ${halFile} ${queryBedFile} ${queryGenome} ${targetGenome} ${outFile}