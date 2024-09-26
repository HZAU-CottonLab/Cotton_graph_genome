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
    # * Turn hal obtained by minigraph into vg.
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
        #* reference J85
        referencePrefix="J85",
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
    # * Extract the coordinates of the breakpoint where PAV occurs.
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
    # TODO Get the real coordinates of SVs breakpoint in each PATH.
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
        cut -f3 {input.bpFile} |grep "^[<>]" |sed 's/[<>]/\\n/g'|sort -n |uniq |sed '/^$/d' >{output.bp_nodeList}
        #* Analyze the real coordinates of node in each path.
        for path in $(vg paths -x {input.pgFile} -L); do
            vg find -x {input.pgFile}  -N  {output.bp_nodeList} -P ${{path}}| \
                awk -v pathname=${{path}} 'NF==2{{print $0"\t"pathname
                    }}NF==1{{print $1"\t-1\t"pathname}}NF>2{{print $1"\t-1\t"pathname}}' >>{output.node_map_postion}
        done
        """


rule bpSite_update_05:
    # * Update the whole SVs file according to the real coordinates of the breakpoint.
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
