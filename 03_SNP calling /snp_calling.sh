# date: 2022-11-01
###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2022-11-04 17:37:40
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 09:50:28
# @@param:
###

module load sentieon/202112
module load SAMtools/1.9

# input

# mkdir SNP_call
sampleList=$1
basedir='/public/home/zpliu/Pan-genome/Genotype_diploid/'
rawBase='/data/cotton/MaojunWang/Raw_Data/DiploidCotton_Resequencing/rawdata/'

#for i in `grep Gb -v  ${sampleList}`; do
for i in $(cat $sampleList); do
    fq1=$rawBase/${i}_R1.fq.gz
    fq2=$rawBase/${i}_R2.fq.gz

    fasta=/public/home/zpliu/Pan-genome/Genotype_diploid/bwa_index/J85.fa

    # output dir
    workdir=$basedir/SNP_call/$i
    [ -d $workdir ] && rm $workdir/*

    # thread
    nt=10

    platform="ILLUMINA"
    mq=30
    [ ! -d $workdir ] && mkdir -p $workdir
    cd $workdir

    rawCram=$i.cram
    sortedCram=$i.q$mq.sorted.cram
    depCram=$i.deduped.cram
    realnCram=$i.realn.cram
    outvcf=$i.vcf

    # queue
    bsub -J ${i} -o %J.${i}.out -e %J.${i}.err -n $nt -q normal "

# *******************************calling ***********************************
exec > $workdir/$i.callVCF.log 2>&1

# 1. Map reads to reference

(sentieon bwa mem -M -R '@RG\tID:${i}\tSM:${i}\tPL:$platform' \
-t $nt -K 10000000 $fasta $fq1 $fq2 || echo -n 'error' )| \
samtools sort -@ $nt --output-fmt CRAM \
--reference $fasta -o $rawCram - && samtools index -@ $nt $rawCram

samtools view -hCS -T $fasta -q $mq -o $sortedCram $rawCram && \
samtools index -@ $nt $sortedCram

samtools flagstat $rawCram > $i.stat.raw.txt && \
samtools flagstat $sortedCram > $i.stat.q$mq.txt &

# 2. Calculate data metrics

sentieon driver -r $fasta -t $nt -i $sortedCram \
--algo MeanQualityByCycle ${i}_mq_metrics.txt \
--algo QualDistribution ${i}_qd_metrics.txt \
--algo GCBias \
--summary ${i}_gc_summary.txt ${i}_gc_metrics.txt \
--algo AlignmentStat \
--adapter_seq '' ${i}_aln_metrics.txt \
--algo InsertSizeMetricAlgo ${i}_is_metrics.txt

sentieon plot metrics -o ${i}_metrics-report.pdf gc=${i}_gc_metrics.txt \
qd=${i}_qd_metrics.txt mq=${i}_mq_metrics.txt isize=${i}_is_metrics.txt

sentieon driver -r $fasta -t $nt -i $sortedCram \
--algo LocusCollector --fun score_info ${i}_score.txt

# 3. Remove duplicates

sentieon driver -r $fasta -t $nt -i $sortedCram --algo Dedup --rmdup \
--cram_write_options version=3.0 \
--score_info ${i}_score.txt --metrics ${i}_dedup_metrics.txt $depCram && \
rm $sortedCram

# 4. Variant calling

sentieon driver -r $fasta -t $nt -i $depCram \
--algo Haplotyper --genotype_model multinomial \
--emit_mode gvcf ${i}.gvcf && \
sentieon util vcfconvert ${i}.gvcf ${i}.gvcf.gz && \
rm ${i}.gvcf && rm $rawCram && touch succeed.txt
"
done
