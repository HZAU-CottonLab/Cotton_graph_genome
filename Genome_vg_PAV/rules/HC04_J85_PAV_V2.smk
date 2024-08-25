Agenome = [
    "A01",
    "A02",
    "A03",
    "A04",
    "A05",
    "A06",
    "A07",
    "A08",
    "A09",
    "A10",
    "A11",
    "A12",
    "A13",
]

rule all:
    input:
        VG_chr_rule02=expand("halMutation/HC04_At_J85/{Chr}_filter.vg", Chr=Agenome),
        VG_deconstruct_rule03=expand(
            "VG_graph/deconstruct_VCF/HC04_At_vs_J85/{Chr}.vcf", Chr=Agenome
        ),
        PAV_data_rule04=expand(
            "VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-info.txt", Chrom=Agenome
        ),
        node_position_rule05=expand(
            "VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_bp_map.txt",Chrom=Agenome
        ),
        breakPoint_rule06=expand(
            "VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-position.txt",Chrom=Agenome
        )


rule A_genome_vg_01:
    input:
        hal="/public/home/zpliu/Pan-genome/Cactus-Pan/subGenomes/Gossypium-steps-output_V2/evolverGossypium.hal",
    output:
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/HC04_At_J85/A_subgenome.pg",
    resources:
        mem_mb=20000,
    threads: 1
    shell:
        """
        module load cactus/2.5.1
        #* 提取A基因组对应树的graph
        hal2vg {input.hal} --chop 32 --rootGenome Anc3 --progress --noAncestors >{output.pgFile}
        """


rule split_VG_by_chroms_02:
    # TODO 按照染色体对应关系，将包含多个物种染色体的graph进行分割
    input:
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/HC04_At_J85/A_subgenome.pg",
        path_select="halMutation/HC04_At_J85/{Chr}_path.txt",
    output:
        Chr_pg="halMutation/HC04_At_J85/{Chr}_filter.vg",
    threads: 2
    resources:
        mem_mb=150000,
    shell:
        """
        module load vg/1.46.0
        vg paths -p {input.path_select} -x  {input.pgFile} -r|vg mod - -N -t {threads} >{output.Chr_pg}
        """


rule vg_deconstruct_02:
    # TODO 解构VG成VCF这一步消耗内存比较大
    input:
        pgFile="halMutation/HC04_At_J85/{Chr}_filter.vg",
    output:
        vcf="VG_graph/deconstruct_VCF/HC04_At_vs_J85/{Chr}.vcf",
    params:
        #* 以J85为作为reference解构graph
        referencePrefix="J85",
        #* 值越大得到的变异越多?
        context_jaccard=50000,
    threads: 4
    resources:
        mem_mb=100000,
    shell:
        """
        module load vg/1.46.0 
        vg deconstruct -e -v -a -c  {params.context_jaccard} -P {params.referencePrefix} {input.pgFile} -t {threads}  >{output.vcf}
        """


rule VCF2PAV_03:
    # * 提取发生PAV的断点坐标, 断点是对应一个graph中的node, 对于该node内的序列仍旧是保守的
    input:
        rawVCF="VG_graph/deconstruct_VCF/HC04_At_vs_J85/{Chrom}.vcf",
    output:
        PAVinfo="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-info.txt",
    params:
        #* 这里的wildcards是A01，需要使用J85_Chr01进行表示reference Path
        referencePathName=lambda wildcards: "J85#J85_Chr{}".format(
            wildcards.Chrom[1:], wildcards.Chrom[1:]
        ),
        VCFskipRows=11,
    threads: 1
    resources:
        mem_mb=10000,
    shell:
        """
        python script/filter_VCF.py {params.VCFskipRows} {input.rawVCF} {params.referencePathName} {output.PAVinfo}
        """

rule Node_in_realPosition_04:
    # TODO 获取SVs断点在每个PATH中的真实坐标
    input:
        bpFile="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-info.txt",
        pgFile="halMutation/HC04_At_J85/{Chrom}_filter.vg",
    output:
        bp_nodeList=temp("VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point_temp.txt"),
        node_map_postion="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_bp_map.txt",
    threads: 1
    resources:
        mem_mb=15000,
    shell:
        """
        module load vg/1.46.0
        #! 这需要修改! 提取断点序列信息
        cut -f3 {input.bpFile} |grep "^[<>]" |sed 's/[<>]/\\n/g'|sort -n |uniq |sed '/^$/d' >{output.bp_nodeList}
        #* 分析node在每个path中的真实坐标
        for path in $(vg paths -x {input.pgFile} -L); do
            vg find -x {input.pgFile}  -N  {output.bp_nodeList} -P ${{path}}| \
                awk -v pathname=${{path}} 'NF==2{{print $0"\t"pathname
                    }}NF==1{{print $1"\t-1\t"pathname}}NF>2{{print $1"\t-1\t"pathname}}' >>{output.node_map_postion}
        done
        """

rule bpSite_update_05:
    # * 根据断点的真实坐标,更新整个SVs文件
    input:
        bp_map="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_bp_map.txt",
        SV_info="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-info.txt",
    output:
        SV_position="VG_graph/PAV_V2/HC04_At_vs_J85/{Chrom}_break_point-position.txt",
    resources:
        mem_mb=15000,
    threads: 1
    shell:
        """
        python script/breakpointMap.py {input.SV_info} {input.bp_map} {output.SV_position}
        """