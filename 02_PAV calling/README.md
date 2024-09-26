<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2024-09-26 10:07:46
 * @LastEditors: zpliu
 * @LastEditTime: 2024-09-26 19:22:29
 * @@param: 
-->

### Cacuts between AD1 and A2

```bash
cactus-prepare seqFile.config --defaultCores 10 \
    --outDir Gossypium-steps-output \
    --outSeqFile Gossypium-steps-output/evolverGossypium.txt \
    --outHal Gossypium-steps-output/evolverGossypium.hal --jobStore jobstore >Gossypium-step-by-step.sh
```

### miniGraph-Cactus in 15 A2 and 35 AD1

```bash 

module load cactus/2.5.1
for chr in {01..13}; do
    mkdir -p panGenome-D${chr}
    for i in $(ls /public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/D_genome/genomeSequence/D${chr}); do
        echo $(echo ${i} | sed 's/\.fasta//g') /public/home/zpliu/Pan-genome/Cactus-Pan/AD_genomes/D_genome/genomeSequence/D${chr}/${i} >>panGenome-D${chr}/genomePangenome.txt
    done
done

#* Chrom
for chr in {01..13}; do
    bsub -q normal -n 36 -e D${chr}.err -o D${chr}.out -J D${chr} "
    module load cactus/2.5.1
    cactus-pangenome jobStore_D${chr} panGenome-D${chr}/genomePangenome.txt --outDir panGenome-D${chr} \
        --outName D${chr} --reference HC04-D${chr} --maxCores 36
    "
done
```


### Directory `graph_vg`
+ `microRearrangement` Get Micro-rearrangement from graph genome
+ `rules` and `script`  Deconstruct PAV from VG file


### Directory `halLiftover` 

+ `bed_liftover.py` Unify the coordinates of multiple genomes



###  Directory `long-readCalling` 

> Scripts for identification of PAV in tetraploid cotton genomes based on ONT long-reads.

+ `01.sniffles_cuteSV_svim.sh`  Three bioinformatics tools, sniffles, cuteSV, and svim, were used to detect SVs based on ONT long-read mapping of 25 AD1 accessions
+ `02.SV merge genotype.sh` SURVIVOR merges SVs from multiple tools and SVs from multiple materials to obtain the final SV Genotype