SamplListFile=config["input"]
File=open(SamplListFile,'r')

Samplelist_stage={
    "4DPA":[],
    "8DPA":[],
    "12DPA":[],
    "16DPA":[],
    "20DPA":[],
}
for line in File:
    sample,stage = line.strip("\n").split("\t")
    Samplelist_stage[stage].append(
        sample
    )
File.close()

rule all:
    input: 
        junctionFile_4DPA=expand( 
            "4DPA/{sample}/gene_expression.txt",sample=Samplelist_stage['4DPA']
        ),
        junctionFile_8DPA=expand( 
            "8DPA/{sample}/gene_expression.txt",sample=Samplelist_stage['8DPA']
        ),
        junctionFile_12DPA=expand( 
            "12DPA/{sample}/gene_expression.txt",sample=Samplelist_stage['12DPA']
        ),
        junctionFile_16DPA=expand( 
            "16DPA/{sample}/gene_expression.txt",sample=Samplelist_stage['16DPA']
        ),
        junctionFile_20DPA=expand( 
            "20DPA/{sample}/gene_expression.txt",sample=Samplelist_stage['20DPA']
        )


rule Expression_gene_4DPA:
    #TODO 从bam文件中提取exon-exon的junction reads
    input: 
        bamFile='/cotton/Liuzhenping/parallel_sQTLs/A2_RNA_bams/4DPA/{sample}_sort_q25.bam',
        gtf='reference.gff3'
    output: 
        expressionFile='4DPA/{sample}/gene_expression.txt'
    threads: 1
    params:
        Ballgowndir='4DPA/{sample}/'
    resources:
        mem_mb=10000
    shell: 
        """
            module load StringTie/2.1.4
            stringtie  {input.bamFile} --rf  \
                -G {input.gtf} -e \
                -p {threads} \
                -A {output.expressionFile} \
                -b {params.Ballgowndir}
        """

rule Expression_gene_8DPA:
    #TODO 从bam文件中提取exon-exon的junction reads
    input: 
        bamFile='/cotton/Liuzhenping/parallel_sQTLs/A2_RNA_bams/8DPA/{sample}_sort_q25.bam',
        gtf='reference.gff3'
    output: 
        expressionFile='8DPA/{sample}/gene_expression.txt'
    threads: 1
    params:
        Ballgowndir='8DPA/{sample}/'
    resources:
        mem_mb=10000
    shell: 
        """
            module load StringTie/2.1.4
            stringtie  {input.bamFile} --rf  \
                -G {input.gtf} -e \
                -p {threads} \
                -A {output.expressionFile} \
                -b {params.Ballgowndir}
        """

rule Expression_gene_12DPA:
    #TODO 从bam文件中提取exon-exon的junction reads
    input: 
        bamFile='/cotton/Liuzhenping/parallel_sQTLs/A2_RNA_bams/12DPA/{sample}_sort_q25.bam',
        gtf='reference.gff3'
    output: 
        expressionFile='12DPA/{sample}/gene_expression.txt'
    threads: 1
    params:
        Ballgowndir='12DPA/{sample}/'
    resources:
        mem_mb=10000
    shell: 
        """
            module load StringTie/2.1.4
            stringtie  {input.bamFile} --rf  \
                -G {input.gtf} -e \
                -p {threads} \
                -A {output.expressionFile} \
                -b {params.Ballgowndir}
        """

rule Expression_gene_16DPA:
    #TODO 从bam文件中提取exon-exon的junction reads
    input: 
        bamFile='/cotton/Liuzhenping/parallel_sQTLs/A2_RNA_bams/16DPA/{sample}_sort_q25.bam',
        gtf='reference.gff3'
    output: 
        expressionFile='16DPA/{sample}/gene_expression.txt'
    threads: 1
    params:
        Ballgowndir='16DPA/{sample}/'
    resources:
        mem_mb=10000
    shell: 
        """
            module load StringTie/2.1.4
            stringtie  {input.bamFile} --rf  \
                -G {input.gtf} -e \
                -p {threads} \
                -A {output.expressionFile} \
                -b {params.Ballgowndir}
        """

rule Expression_gene_20DPA:
    #TODO 从bam文件中提取exon-exon的junction reads
    input: 
        bamFile='/cotton/Liuzhenping/parallel_sQTLs/A2_RNA_bams/20DPA/{sample}_sort_q25.bam',
        gtf='reference.gff3'
    output: 
        expressionFile='20DPA/{sample}/gene_expression.txt'
    threads: 1
    params:
        Ballgowndir='20DPA/{sample}/'
    resources:
        mem_mb=10000
    shell: 
        """
            module load StringTie/2.1.4
            stringtie  {input.bamFile} --rf  \
                -G {input.gtf} -e \
                -p {threads} \
                -A {output.expressionFile} \
                -b {params.Ballgowndir}
        """
# rule Intron_cluster:
    #TODO 对intron进行聚类

