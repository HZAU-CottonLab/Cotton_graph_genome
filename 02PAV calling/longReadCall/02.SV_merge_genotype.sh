=====================SV merge=======================
#1.将三个软件的结果取两个结果一致的，去掉三个结果不一样的。
module load SURVIVOR/1.0.6 

for line in `ls ../svim/*vcf`; 
do name=`basename $line| sed s/.vcf//`;
printf "$name.cuteSV.INS.vcf\n$name.svim.INS.vcf\n$name.sniffles.INS.vcf\n" > $name.INS.txt;
SURVIVOR  merge $name.INS.txt 1000 2 1 -1 -1 50 merged_$name.INS.vcf;
python consenus.py merged_$name.INS.vcf
done;

for line in `ls ../svim/*vcf`; 
do name=`basename $line| sed s/.vcf//`;
printf "$name.cuteSV.DEL.vcf\n$name.svim.DEL.vcf\n$name.sniffles.DEL.vcf\n" > $name.DEL.txt;
SURVIVOR  merge $name.DEL.txt 1000 2 1 -1 -1 50 merged_$name.DEL.vcf;
python consenus.py merged_$name.DEL.vcf
done;

lh consenus_*DEL* |awk '{print $9}' > DEL.consenus.txt 
lh consenus_*INS* |awk '{print $9}' > INS.consenus.txt

#2. 合并所有材料的变异位点
SURVIVOR merge DEL.consenus.txt 1000 1 1 -1 -1 -1 ../genotype/DEL.vcf
SURVIVOR merge INS.consenus.txt 1000 1 1 -1 -1 -1 ../genotype/INS.vcf

=====================SV Genotype=======================
grep '##' DEL.vcf > DEL.site.vcf
grep -v "##" DEL.vcf |cut -f 1-9 >> DEL.site.vcf

grep '##' INS.vcf > INS.site.vcf
grep -v '##' INS.vcf |cut -f 1-9 >> INS.site.vcf

module load sniffles/2.0.2

#sniffles --input ../DC001_reads_mapping_ref.sorted.bam -t 4 --vcf DC001.DEL.vcf --genotype-vcf DEL.site.vcf

for line in ../*bam;do
name=`basename $line | sed 's/_reads_mapping_ref.sorted.bam//g'`
bsub -n 4 -q q2680v2 " sniffles --input $line -t 4 --vcf ${name}.INS.vcf --genotype-vcf INS.site.vcf  "
bsub -n 4 -q q2680v2 " sniffles --input $line -t 4 --vcf ${name}.DEL.vcf --genotype-vcf DEL.site.vcf  "
done


lh DC*DEL.vcf |awk '{print $9}' > DEL.lst.txt
module load SURVIVOR/1.0.6
SURVIVOR merge DEL.lst.txt 1000 1 1 -1 -1 1 final.DEL.vcf

lh DC*INS.vcf |awk '{print $9}' > INS.lst.txt
SURVIVOR merge INS.lst.txt  1000 1 1 -1 -1 1 final.INS.vcf
