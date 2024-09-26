#bwa index
bwa index HC04.ctg.fasta
python generate_site_positions.py DpnII HC04 HC04.ctg.fasta
awk 'BEGIN{OFS="\t"}{print $1, $NF}' P20_DpnII.txt > HC04.chrom.sizes

#juicer
bash /public/home/jyli/software/juicer/scripts/juicer.sh -g HC04 -d /public/home/jyli/LJY_GhONT_Project/HiFi_AS/HC04_3Cell_HiC/juicer/ 
	-D /public/home/jyli/software/juicer/CPU -z /public/home/jyli/LJY_GhONT_Project/HiFi_AS/HC04_3Cell_HiC/juicer/references/HC04.ctg.fasta 
	-y /public/home/jyli/LJY_GhONT_Project/HiFi_AS/HC04_3Cell_HiC/juicer/references/HC04_DpnII.txt 
	-p /public/home/jyli/LJY_GhONT_Project/HiFi_AS/HC04_3Cell_HiC/juicer/references/HC04.chrom.sizes -t 30 
	
#3d-dna
run-asm-pipeline.sh -m haploid -r 0 HC04.ctg.fasta /data/cotton/TMP/LiJianYing_ONTs/HiFi_AS/Juicer/HC04/aligned/merged_nodups.txt
run-asm-pipeline-post-review.sh -r HC04.ctg.0.review.assembly HC04.ctg.fasta /data/cotton/TMP/LiJianYing_ONTs/HiFi_AS/Juicer/HC04/aligned/merged_nodups.txt

#RagTag
for ID in HC15 HW03 HW05 HW06 HW07 P01 P02 P04 P19 P20 TW007 TW013 TW026 TW029 TW031 TW055 TW064 TW075 TW077 TW091 TW094 TW100 TW134 XJ74 XZ142 ZY006 ZY10 ZY184 ZY236 ZY238 ZY354 ZY381 ZY384 ZY461
do
ragtag.py scaffold HC04.Chr.fa ${ID}.asm.bp.p_ctg.fasta -t 10 -o ${ID}
done

