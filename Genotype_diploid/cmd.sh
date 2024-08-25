###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-02-17 09:00:32
 # @LastEditors: zpliu
 # @LastEditTime: 2023-03-10 14:51:31
# @@param:
###
#* 复制参考基因组序列
cp /cotton/JianyingLi/PAN_Graph/A2/01.HiFi_mapping_J85/J85.fa .
bsub -q normal -n 1 -J bwa -e bwa_index.err -o bwa_index.out "
    bwa index J85.fa
"

#* 合并所有的gvcf文件
module load sentieon/202112
bsub -q normal -n 10 -e test.err -o test.out -J sentieon "
    cat gvcf_list.txt|sentieon driver -t 10 -r /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa --algo GVCFtyper J85-14Sample.vcf.gz -
"

#* 将gvcf转为vcf
module load Singularity/3.7.3
inputPath='/public/home/zpliu/Pan-genome/Genotype_diploid/SNP_call/'
for sample in $(cat All_samplelist.txt Sample_with_bam.txt | sort | uniq -u); do
    bsub -q q2680v2 -n 4 -e ${sample}.err -o ${sample}.out "
        singularity exec --nv $IMAGE/clara-parabricks/4.0.1-1.sif   pbrun genotypegvcf \
                    --ref /public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa \
                    --in-gvcf ${inputPath}/${sample}/${sample}.g.vcf.gz \
                    --out-vcf ${inputPath}/${sample}/${sample}.vcf \
                    --num-threads 4 
"
done
#TODO 压缩建索引
inputPath='/public/home/zpliu/Pan-genome/Genotype_diploid/SNP_call/'
for sample in $(cat All_samplelist.txt); do
    bsub -q q2680v2 -n 1 -e ${sample}.err -o ${sample}.out "
        bgzip -c -f ${inputPath}/${sample}/${sample}.vcf
        bcftools index ${inputPath}/${sample}/${sample}.vcf.gz 
    "
done



#TODO 合并所有的vcf文件
