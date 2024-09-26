###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-02-27 16:29:37
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 09:49:31
 # @@param: 
### 
#* Removing PCR duplication and indexing 
module load picard/2.23.9
outPath=$1
sampleName=$2
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
        INPUT=${outPath}/${sampleName}_sorted_q20.bam \
        OUTPUT=${outPath}/${sampleName}_srt_q20_redup.bam \
        METRICS_FILE=${outPath}/${sampleName}_metrics.txt 

#* Build index
java -jar ${EBROOTPICARD}/picard.jar BuildBamIndex  INPUT=${outPath}/${sampleName}_srt_q20_redup.bam