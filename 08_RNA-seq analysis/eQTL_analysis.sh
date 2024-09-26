###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2024-09-26 11:32:19
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 11:32:20
 # @@param: 
### 
for stage in  4DPA 12DPA 16DPA 20DPA 8DPA; do
        plink_prefix='HC04_snp_filter_maf'
        softPath='/public/home/zpliu/LZP_sQTL_A2_At/leafCutter/AD1_HC04/script/'
        phenotype="${stage}_fastqtl.bed.gz"
        covariant="${stage}_covariant.txt"
        for i in $(ls splitList | grep split_$stage); do
                bsub -q normal -n 1 -e ${stage}.err -o ${stage}.out -J AD1_${i} -M 10G "
        python ${softPath}/QTL_mapping.py ${plink_prefix} ${phenotype} ${covariant} ${stage}_${i}_trans.txt t splitList/${i}
        "
        done
done