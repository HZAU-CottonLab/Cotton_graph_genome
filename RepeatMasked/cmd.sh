###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-04-27 20:32:12
 # @LastEditors: zpliu
 # @LastEditTime: 2023-05-30 21:31:57
# @@param:
###

#TODO 运行snakemake流程
#* 测试一下流程是否能够完整跑通
module load RepeatMasker/4.1.0

snakemake  --configfile inputConfig.json \
   -np -s mask.smk

#TODO 在集群中运行流程
snakemake --profile lsf \
    --configfile inputConfig.json \
    -j 50 -s mask.smk

#TODO 使用EDTA hard-masked repeat后的基因组序列进行TRF的masked 
#! masked 采用的是hard-masked
#? 最终masked后的文件：  Finally_Masked_HC04.fa


#TODO 进行segment duplicate注释, 这个得手动交任务
#? 结果文件
HC04_SD.annotation.sd
HC04_SD.annotation.sd.elem.txt 


#*seqtk是0-base的，坐标需要减1
cut -f1-3 Masked/HC04_trf.bed |bedtools sort -i - | bedtools merge -i - \
    |awk '{print $1"\t"$2-1"\t"$3}' \
    |/cotton/Liuzhenping/Pan-genome/software/seqtk-1.3/seqtk  seq -l 60 -M /dev/stdin HC04-softMasked.fa >HC04-Repeat_TRF_softMasked.fa


#TODO 获取EDTA和TRF注释的区间
python ./scripts/get_softMaskedBed.py HC04-Repeat_TRF_softMasked.fa HC04_maskedRegion.bed 


#! 将初步过滤后的SD序列与Repeat区域取交集
cat HC04_SD_filter.txt|cut -f1-3,11 \
    |intersectBed  -a - -b HC04_Repeat_TRF_maskedRegion.bed   -f 0.7  -wa -wb  >intersectData.txt

cat HC04_SD_filter.txt|cut -f4-6,11 \
    |intersectBed  -a - -b HC04_Repeat_TRF_maskedRegion.bed   -f 0.7 -wb -wa   >intersectData_V2.txt

#! 根据筛选结果从SD文件中保留的最终的SD集合





#TODO 使用许忠平构建的J668的TE lib序列进行Repeatmask
snakemake --configfile inputConfig.json \
    -s mask.smk -np --cores 1 --target  RepeatSeq_RepeatMasker_soft_01:


