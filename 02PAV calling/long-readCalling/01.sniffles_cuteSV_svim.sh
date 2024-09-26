#minimap2 + sniffles
minimap2 -t 12 -R '@RG\tID:DC113\tSM:DC113' -ax map-hifi J85.fa reads/DC113.fq.gz | sentieon util sort -r J85.fa --sam2bam -t  12 -o DC113.bam -i - && sniffles --input DC113.bam --snf DC113.SNF -t 12 --vcf DC113.vcf

#cuteSV
for line in ../*bam;do
name=`basename $line | sed s/_reads_mapping_ref.sorted.bam//`
mkdir $name
bsub -n 12 -e $name.err -o $name.out  "/public/home/xhe/miniconda3/envs/nanoSV/bin/cuteSV $line ../J85.fa $name.vcf -t 12 $name \
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
--max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 \
--sample $name --min_size 50 --genotype "
done

#svim
for line in ../*bam;do
name=`basename $line | sed s/_reads_mapping_ref.sorted.bam//`;
mkdir $name
bsub -n 1 -e $name.err -o $name.out  "/public/home/xhe/miniconda3/envs/nanoSV/bin/svim alignment $name $line ../J85.fa "
echo $name;
done