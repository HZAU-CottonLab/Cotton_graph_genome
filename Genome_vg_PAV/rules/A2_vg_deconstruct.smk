J85Chroms = [
    "Chr01",
    "Chr02",
    "Chr03",
    "Chr04",
    "Chr05",
    "Chr06",
    "Chr07",
    "Chr08",
    "Chr09",
    "Chr10",
    "Chr11",
    "Chr12",
    "Chr13",
]


rule all:
    input:
        pgFile_rule01=expand("halMutation/J85_vs_A2s/{Chrom}.pg", Chrom=J85Chroms),
        vcf_rule02=expand(
            "VG_graph/deconstruct_VCF/J85_vs_A2s/{Chrom}.vcf", Chrom=J85Chroms
        ),
        PAVinfo_rule03=expand(
            "VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-info.txt", Chrom=J85Chroms
        ),
        node_map_postion_rule04=expand(
            "VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_bp_map.txt", Chrom=J85Chroms
        ),
        SV_position_rule05=expand(
            "VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-position.txt",
            Chrom=J85Chroms,
        ),


rule hal2vg_J85_01:
    # * 将minigraph得到的hal文件转为vg的graph文件
    input:
        halFile="/public/home/zpliu/Pan-genome/Cactus-Pan/A2_genome/panGenome-{Chrom}/{Chrom}.full.hal",
    output:
        pgFile="halMutation/J85_vs_A2s/{Chrom}.pg",
    threads: 1
    resources:
        mem_mb=100000,
    shell:
        """
        module load cactus/2.5.1
        hal2vg {input.halFile} --chop 32 --rootGenome Anc0 --progress \
            --noAncestors --ignoreGenomes _MINIGRAPH_ >{output.pgFile}
        """


rule J85_vg_deconstruct_02:
    input:
        pgFile="halMutation/J85_vs_A2s/{Chrom}.pg",
    output:
        vcf="VG_graph/deconstruct_VCF/J85_vs_A2s/{Chrom}.vcf",
    params:
        #* 以J85为作为reference解构graph
        referencePrefix="J85",
        #* 值越大得到的变异越多?
        context_jaccard=50000,
    threads: 4
    resources:
        mem_mb=10000,
    shell:
        """
        module load vg/1.46.0 
        vg deconstruct -e -v -a -c  {params.context_jaccard} -P {params.referencePrefix} {input.pgFile} -t {threads}  >{output.vcf}
        """


rule J85_VCF2PAV_03:
    # * 提取发生PAV的断点坐标, 断点是对应一个graph中的node, 对于该node内的序列仍旧是保守的
    input:
        rawVCF="VG_graph/deconstruct_VCF/J85_vs_A2s/{Chrom}.vcf",
    output:
        PAVinfo="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-info.txt",
    params:
        referencePathName=lambda wildcards: "J85-{}#{}".format(
            wildcards.Chrom, wildcards.Chrom
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
        bpFile="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-info.txt",
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/J85_vs_A2s/{Chrom}.pg",
    output:
        bp_nodeList=temp("VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point_temp.txt"),
        node_map_postion="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_bp_map.txt",
    threads: 1
    resources:
        mem_mb=10000,
    shell:
        """
        module load vg/1.46.0
        #! 这需要修改
        cut -f3 {input.bpFile} |grep "^[<>]" |sed 's/[<>]/\\n/g'|sort -n |uniq |sed '/^$/d' >{output.bp_nodeList}
        #* 分析node在每个path中的真实坐标
        for path in $(vg paths -x {input.pgFile} -L); do
            vg find -x {input.pgFile}  -N  {output.bp_nodeList} -P ${{path}}| \
                awk -v pathname=${{path}} 'NF==2{{print $0"\t"pathname
                    }}NF==1{{print $1"\t-1\t"pathname}}NF>2{{print $1"\t-1\t"pathname}}' >>{output.node_map_postion}
        done
        """


rule bpSite_update_05:
    # * 根据断点的真实坐标，更新整个SVs文件
    input:
        bp_map="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_bp_map.txt",
        SV_info="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-info.txt",
    output:
        SV_position="VG_graph/PAV_V2/J85_vs_A2s/{Chrom}_break_point-position.txt",
    resources:
        mem_mb=10000,
    threads: 1
    shell:
        """
        python script/breakpointMap.py {input.SV_info} {input.bp_map} {output.SV_position}
        """
