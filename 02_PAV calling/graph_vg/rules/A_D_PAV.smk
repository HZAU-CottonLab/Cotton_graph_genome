Agenome = [
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
        VG_chr_rule02=expand(
            "/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/{Chr}_filter.vg",
            Chr=Agenome,
        ),
        vcf_rule03=expand(
            "VG_graph/deconstruct_VCF/A_vs_D/{Chr}.vcf",
            Chr=Agenome,
        ),
        PAVinfo_rule04=expand(
            "VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-info.txt",
            Chr=Agenome,
        ),
        node_map_postion_rule05=expand(
            "VG_graph/PAV_V2/A_vs_D/{Chr}_bp_map.txt", Chr=Agenome
        ),
        SV_position_rule05=expand(
            "VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-position.txt",
            Chr=Agenome,
        ),


rule hal_2_vg:
    input:
        hal="/public/home/zpliu/Pan-genome/Cactus-Pan/subGenomes/Gossypium-steps-output_V2/evolverGossypium.hal",
    output:
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/A_D_subgenome.pg",
    resources:
        mem_mb=250000,
    threads: 1
    shell:
        """
        module load cactus/2.5.1
        hal2vg {input.hal} --chop 32 --rootGenome Anc1 --progress \
            --ignoreGenomes A1,A1a \
            --noAncestors >{output.pgFile}
        """


rule split_VG_by_chroms_02:
    input:
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/A_D_subgenome.pg",
        path_select="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/{Chr}_path.txt",
    output:
        Chr_pg="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/{Chr}_filter.vg",
    threads: 2
    resources:
        mem_mb=200000,
    shell:
        """
        module load vg/1.46.0
        vg paths -p {input.path_select} -x  {input.pgFile} -r|vg mod - -N -t {threads} >{output.Chr_pg}
        """


rule vg_deconstruct_03:
    input:
        pgFile="/public/home/zpliu/Pan-genome/SV_parallele_V2/halMutation/A_vs_D/{Chr}_filter.vg",
    output:
        vcf="VG_graph/deconstruct_VCF/A_vs_D/{Chr}.vcf",
    params:
        referencePrefix="J85",
        context_jaccard=50000,
    threads: 4
    resources:
        mem_mb=100000,
    shell:
        """
        module load vg/1.46.0 
        vg deconstruct -e -v -a -c  {params.context_jaccard} -P {params.referencePrefix} {input.pgFile} -t {threads}  >{output.vcf}
        """


rule VCF2PAV_04:
    input:
        rawVCF="VG_graph/deconstruct_VCF/A_vs_D/{Chr}.vcf",
    output:
        PAVinfo="VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-info.txt",
    params:
        referencePathName=lambda wildcards: "J85#J85_{}".format(wildcards.Chr),
        VCFskipRows=11,
    threads: 1
    resources:
        mem_mb=15000,
    shell:
        """
        python script/filter_VCF.py {params.VCFskipRows} {input.rawVCF} {params.referencePathName} {output.PAVinfo}
        """


rule Node_in_realPosition_05:
    input:
        bpFile="VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-info.txt",
        pgFile="halMutation/A_vs_D/{Chr}_filter.vg",
    output:
        bp_nodeList=temp("VG_graph/PAV_V2/A_vs_D/{Chr}_break_point_temp.txt"),
        node_map_postion="VG_graph/PAV_V2/A_vs_D/{Chr}_bp_map.txt",
    threads: 1
    resources:
        mem_mb=15000,
    shell:
        """
        module load vg/1.46.0
        cut -f3 {input.bpFile} |grep "^[<>]" |sed 's/[<>]/\\n/g'|sort -n |uniq |sed '/^$/d' >{output.bp_nodeList}
        for path in $(vg paths -x {input.pgFile} -L); do
            vg find -x {input.pgFile}  -N  {output.bp_nodeList} -P ${{path}}| \
                awk -v pathname=${{path}} 'NF==2{{print $0"\t"pathname
                    }}NF==1{{print $1"\t-1\t"pathname}}NF>2{{print $1"\t-1\t"pathname}}' >>{output.node_map_postion}
        done
        """


rule bpSite_update_05:
    input:
        bp_map="VG_graph/PAV_V2/A_vs_D/{Chr}_bp_map.txt",
        SV_info="VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-info.txt",
    output:
        SV_position="VG_graph/PAV_V2/A_vs_D/{Chr}_break_point-position.txt",
    resources:
        mem_mb=15000,
    threads: 1
    shell:
        """
        python script/breakpointMap.py {input.SV_info} {input.bp_map} {output.SV_position}
        """
