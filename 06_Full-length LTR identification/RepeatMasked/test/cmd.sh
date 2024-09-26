samtools faidx /cotton/Liuzhenping/Pan-genome/HC04/HC04_V2/HC04_chr_adjust.fa HC04_D08:37398000-39805600  >HC04_D08_CENH3.fa
/cotton/Liuzhenping/Pan-genome/software/trf409.linux64 HC04_D08_CENH3.fa 2 7 7 80 10 50 1000 -l 25 -h -ngs > HC04_D08_CENH3.dat
python ../scripts/trf2bed_V2.py  HC04_D08_CENH3.dat  HC04_D08_CENH3.bed
