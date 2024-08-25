AD1Samples = config["AD1"]
A2Samples = config["A2"]


rule all:
    input:
        At_genotype_rule01=expand(
            "superPan/At_genotypes/genotype/{sample}_PAV.vcf", sample=AD1Samples
        ),
        Dt_genotype_rule01=expand(
            "superPan/Dt_genotypes/genotype/{sample}_PAV.vcf", sample=AD1Samples
        ),
        A2_genotype_rule01=expand(
            "superPan/A2_genotypes/genotype/{sample}_PAV.vcf", sample=A2Samples
        ),


rule SV_genotype_01:
    # TODO 对Cactus鉴定的PAV位点进行genotype
    input:
        SV_location_VCF="superPan/At_genotypes/Cactue_PAV_site.vcf",
        bamFile="/cotton/JianyingLi/PAN_Graph/AD1/01.HiFi_mapping_HC04/{sample}_reads_mapping_ref.sorted.bam",
    output:
        SV_genotype="superPan/At_genotypes/genotype/{sample}_PAV.vcf",
    threads: 2
    resources:
        mem_mb=20000,
    shell:
        """
        #* 使用V2.0.7版本
        sniffles --input {input.bamFile} -t {threads} --vcf {output.SV_genotype} \
            --genotype-vcf {input.SV_location_VCF} 
        """


rule SV_genotype_V2_01:
    # TODO 对Cactus鉴定的PAV位点进行genotype
    input:
        SV_location_VCF="superPan/Dt_genotypes/Cactue_PAV_site.vcf",
        bamFile="/cotton/JianyingLi/PAN_Graph/AD1/01.HiFi_mapping_HC04/{sample}_reads_mapping_ref.sorted.bam",
    output:
        SV_genotype="superPan/Dt_genotypes/genotype/{sample}_PAV.vcf",
    threads: 2
    resources:
        mem_mb=20000,
    shell:
        """
        #* 使用V2.0.7版本
        sniffles --input {input.bamFile} -t {threads} --vcf {output.SV_genotype} \
            --genotype-vcf {input.SV_location_VCF} 
        """


rule SV_genotype_V3_01:
    # TODO 对Cactus鉴定的J85 PAV位点进行genotype
    input:
        SV_location_VCF="superPan/A2_genotypes/Cactue_PAV_site.vcf",
        bamFile="/cotton/JianyingLi/PAN_Graph/A2/01.HiFi_mapping_J85/BAMs/{sample}_reads_mapping_ref.sorted.bam",
    output:
        SV_genotype="superPan/A2_genotypes/genotype/{sample}_PAV.vcf",
    threads: 2
    resources:
        mem_mb=20000,
    shell:
        """
        #* 使用V2.0.7版本
        sniffles --input {input.bamFile} -t {threads} --vcf {output.SV_genotype} \
            --genotype-vcf {input.SV_location_VCF} 
        """
