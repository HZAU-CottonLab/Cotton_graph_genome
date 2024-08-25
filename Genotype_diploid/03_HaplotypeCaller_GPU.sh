###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-03-04 09:48:03
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-09 21:15:04
 # @@param: 
### 

module load Singularity/3.7.3
reference='/public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa'
SampleName=$1
inputPath='/public/home/zpliu/Pan-genome/Genotype_diploid/SNP_call/'

#* 在142节点提交任务
bsub -q gpu  -gpu "num=1:gmem=12G" -m gpu01 -J ${SampleName} -e ${SampleName}.err -o ${SampleName}.out -n 5 "
    singularity exec --nv  $IMAGE/clara-parabricks/4.0.1-1.sif \
        pbrun haplotypecaller  --ref ${reference}  \
        --in-bam ${inputPath}/${SampleName}/${SampleName}_srt_q20_redup.bam \
        --out-variants ${inputPath}/${SampleName}/${SampleName}.vcf 
"


