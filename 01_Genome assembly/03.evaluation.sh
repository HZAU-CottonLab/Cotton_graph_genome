#LAI 
ltr_finder -D 20000 -d 1000 -L 3500 -l 100 -p 30 -C -M 0.8 ${ID}.chr.fa > ${ID}.chr.finder.scn
LTR_retriever -threads 30 -genome ${ID}.chr.fa -infinder ${ID}.chr.finder.scn
LAI -t 30 -genome ${ID}.chr.fa -intact ${ID}.chr.fa.pass.list -all ${ID}.chr.fa.out

#BUSCO
busco -i ${ID}.chr.fa -l embryophyta_odb10 -o HW07_out -m genome -c 10

#Merqury
meryl k=21 count output read.meryl ${ID}.ccs.fastq.gz
merqury.sh read.meryl ${ID}.asm.bp.p_ctg.fasta ${ID}

