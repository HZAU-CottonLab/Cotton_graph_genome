"""
Descripttion: 
version: 
Author: zpliu
Date: 2023-04-27 19:59:56
LastEditors: zpliu
LastEditTime: 2023-04-27 20:07:11
@param: 
"""
import os
import sys
import pysam
import pandas as pd

DEBUG = True  # * 暂时不删除临时文件


def tempd(filename):
    if DEBUG:
        return filename
    return temp(filename)


FASTA = os.path.abspath(config["fasta"])
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. try samtools faidx {FASTA}"

# * 指定样本名
SM = "HC04"
if "sample" in config:
    SM = config["sample"]

THREADS = 16
if "threads" in config:
    THREADS = config["threads"]

ChromosomeList = ["A{:02}".format(i) for i in range(1, 14)] + [
    "D{:02}".format(i) for i in range(1, 14)
]
if "Chrs" in config:
    ChromosomeList = config["Chrs"]

# * 对于多个材料的基因组
SMS = [SM]
# * 最终生成的sd预测文件


rule all:
    input:
        trfBed=f"Masked/{SM}_trf.bed",
        # MaskedFasta=f"Finally_Masked_HC04.fa",
        sdAnno=f"HC04_SD.annotation_V2.sd",


rule split_fasta:
    input:
        fasta=FASTA,
    output:
        fastas=tempd("Masked/temp/{sample}_{Chr}.fa"),
    resources:
        #* 指8Gb？
        mem_mb=10000,
    run:
        # TODO 将基因组按染色体切割
        fasta = pysam.FastaFile(input.fasta)
        outs = [open(f, "w+") for f in output.fastas]
        outindx = 0
        for name in fasta.references:
            seq = fasta.fetch(name)
            outs[outindx].write(">{}\n{}\n".format(name, seq))
            outindx += 1
        # * 文件读取完成
        for out in outs:
            out.close()


# ---------------------------------------------------------------
############################run TRF##############################
# ---------------------------------------------------------------
rule run_trf:
    input:
        fasta="Masked/temp/{sample}_{Chr}.fa",
    output:
        dat="Masked/temp/{sample}_{Chr}.dat",
    resources:
        mem_mb=10000,
    threads: 1
    shell:
        """
        /cotton/Liuzhenping/Pan-genome/software/trf409.linux64 {input.fasta} 2 7 7 80 10 50 15 -l 25 -h -ngs > {output.dat}
        """


rule trf_2_bed:
    input:
        dats=expand("Masked/temp/{sample}_{Chr}.dat", sample=SMS, Chr=ChromosomeList),
    output:
        bed=f"Masked/{SM}_trf.bed",
    script:
        "./scripts/trf2bed.py"


# ---------------------------------------------------------------
############################REPEAT MASKER########################
# TODO 使用EDTA对repeat进行了注释
#! 直接在EDTA masked的基础上进一步对TRF进行masked
# ---------------------------------------------------------------
rule EDTARepeatmasker:
    input:
        fasta=FASTA,
    output:
        output="EDTA/HC04_chr_adjust.fa.mod.MAKER.masked",
    resources:
        mem_mb=30000,
    threads: 36
    shell:
        """
        module load Singularity/3.1.1
        module load Perl/5.26.2 
        cp {input.fasta} $(dirname {output.output})
        cd $(dirname {output.output})
        echo `pwd` $(basename input.fasta)
        singularity exec /public/home/zpliu/TIP/TE_annotion/EDTA.sif EDTA.pl \
            --genome $(basename {input.fasta}) \
            --step all --overwrite 1 --sensitive 0 --anno 1 -t {threads}
        """


#! 直接在EDTA masked的基础上进一步对TRF进行hard-masked
# ? 参考 https://github.com/mrvollger/assembly_workflows/tree/v0.3-Zenodo/workflows
# ? software: https://github.com/lh3/seqtk
rule trf_masked:
    input:
        ReMaskedfasta="EDTA/HC04_chr_adjust.fa.mod.MAKER.masked",
        TRFbed=rules.trf_2_bed.output,
    output:
        finallyMaskedFasta="Finally_Masked_HC04.fa",
    params:
        seqtk="/cotton/Liuzhenping/Pan-genome/software/seqtk-1.3/seqtk",
    threads: 1
    resources:
        mem_mb=10000,
    shell:
        #! 由于seqtk是0-base的起始坐标需要减一,这里还没减
        """
        module load SAMtools/1.9
        cut -f1-3 {input.TRFbed}| bedtools sort -i - | bedtools merge -i - \
            | {params.seqtk}  seq -l 60 -n N -M /dev/stdin {input.ReMaskedfasta} >{output.finallyMaskedFasta}
        samtools faidx {output.finallyMaskedFasta}
        """


# TODO 对masked后的参加基因组进行重复片段的注释
# ? software: https://github.com/0xTCG/biser/


rule segment_annotation:
    input:
        maskedGenome="HC04-Repeat_TRF_softMasked.fa",
    output:
        outPut="HC04_SD.annotation_V2.sd",
        elem="HC04_SD.annotation_V2.sd.elem.txt",
    threads: 20
    resources:
        mem_mb=50000,
    shell:
        """
        python -m biser -o {output.outPut}  \
                --threads {threads}  --gc-heap 8G  --keep-contigs  {input.maskedGenome}
        """


rule RepeatSeq_RepeatMasker_soft_01:
    input:
        DeNovo_Library="/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker/Jin668_allRepeats.lib",
        genome="/public/home/zpliu/Pan-genome/RepeatMasked/HC04-Repeat_TRF_softMasked.fa",
    output:
        soft_tbi="/public/home/zpliu/Pan-genome/RepeatMasked/RepeatMasker/HC04-Repeat_TRF_softMasked.tbi",
    params:
        soft_outdir=lambda w, output: os.path.dirname(output.soft_tbi),
    threads: 20
    log:
        "logs/RepeatSeq_RepeatMasker_soft_01.log",
    resources:
        mem_mb=50000,
    shell:
        """
        module load RepeatModeler/2.0.1
        module load RepeatMasker/4.1.0
        RepeatMasker -pa {threads} -e ncbi -lib {input.DeNovo_Library} \
                -a -s -xsmall -ace -xm -poly -html -gff -u -dir {params.soft_outdir} {input.genome} 2>{log} 1>&2
        """
