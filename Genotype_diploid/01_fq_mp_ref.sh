###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-02-27 14:53:48
 # @LastEditors: zpliu
 # @LastEditTime: 2023-02-27 16:55:21
# @@param:
###


module load SAMtools/1.9
module load BWA/0.7.17

sampleList=$1
basedir='/public/home/zpliu/Pan-genome/Genotype_diploid/'
rawBase='/data/cotton/MaojunWang/Raw_Data/DiploidCotton_Resequencing/rawdata/'
bwaIndex='/public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa'
outBase='/public/home/zpliu/Pan-genome/Genotype_diploid/SNP_call/' 
nthreads=5
platform='ILLUMINA'
for sampleName in $(cat ${sampleList}); do
    outPath=${outBase}/${sampleName}/
    fq1=$rawBase/${sampleName}_R1.fq.gz
    fq2=$rawBase/${sampleName}_R2.fq.gz
    #* 序列比对
    bsub -q normal -n ${nthreads} -J ${sampleName} -e ${sampleName}.err -o ${sampleName}.out "
    bwa mem -R '@RG\tID:${sampleName}\tSM:${sampleName}\tPL:$platform' \
            -t ${nthreads} ${bwaIndex} ${fq1} ${fq2} \
            |samtools  view -bS  -q 20 -@ ${nthreads} \
            |samtools sort -@ ${nthreads}  -o ${outPath}/${sampleName}_sorted_q20.bam 
    "
done
