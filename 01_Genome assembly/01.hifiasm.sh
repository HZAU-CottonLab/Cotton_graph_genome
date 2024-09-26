hifiasm -o ${ID}.asm -t 20 ${ID}.ccs.fastq.gz
awk '{if($1 ~/S/)print ">"$2"\n"$3}' ${ID}.asm.bp.p_ctg.gfa > ${ID}.asm.bp.p_ctg.fasta
seqkit stat ${ID}.asm.bp.p_ctg.fasta -a > ${ID}.asm.bp.p_ctg.fasta.stat