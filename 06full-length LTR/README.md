<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2024-09-26 10:28:35
 * @LastEditors: zpliu
 * @LastEditTime: 2024-09-26 12:14:11
 * @@param: 
-->
### RUN EDTA maksed using Snakemake workflow

```bash
#TODO Run the snakemake process
snakemake --profile lsf \
    --configfile inputConfig.json \
    -j 50 -s mask.smk

#TODO Using the genome sequence after EDTA hard-masked repeat to masked TRF.


#* 0-base
cut -f1-3 Masked/HC04_trf.bed |bedtools sort -i - | bedtools merge -i - \
    |awk '{print $1"\t"$2-1"\t"$3}' \
    |/cotton/Liuzhenping/Pan-genome/software/seqtk-1.3/seqtk  seq -l 60 -M /dev/stdin HC04-softMasked.fa >HC04-Repeat_TRF_softMasked.fa


#TODO get softmask region
python ./scripts/get_softMaskedBed.py HC04-Repeat_TRF_softMasked.fa HC04_maskedRegion.bed 

```



+ `getLTRCluster.py`   Extract the clustering results of full-length LTR from cdhit
+ `run_flLTR_family.py`  Calculate the similarity of LTR between each two genomes
+ `Pan_1000times.py`     Calculate the number of pan-LTRs and core-LTRs based on randomly sampled 1000 times
+ `Mean_SD.py`           Calculate the mean and standard error

