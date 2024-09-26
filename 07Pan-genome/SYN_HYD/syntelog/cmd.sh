###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2023-10-17 17:47:26
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 11:26:02
# @@param:
###
#TODO Ortholog finder
 orthofinder -f {params.inputSeqData_path} -S diamond -M msa -T raxml-ng -t {threads} -a {threads}



#TODO Extract the ID of homologous gene information in the interval.

for i in $(ls D5_split); do
    bsub -q normal -n 1 -e test.err -o test.out -J D5_${i} " 
     python D5_liftover_gene.py D5_split/${i}  D5_split/${i}_out 
    "
done


for i in $(ls At_Dt_split); do
    bsub -q normal -n 1 -e At_Dt_split/test.err -o At_Dt_split/test.out -J A2_D5_${i} " 
     python At_Dt_mosaic_SYN.py At_Dt_split/${i}  At_Dt_split/${i}_out 
    "
done

for i in $(ls At_Dt); do
    bsub -q normal -n 1 -e At_Dt/test.err -o At_Dt/test.out -J At_Dt_${i} " 
     python mosaic_OrthoGeneList.py  At_Dt/${i}  At_Dt/${i}_out 
    "
done


for i in $(ls D5_Dt); do
    bsub -q normal -n 1 -e D5_Dt/test.err -o D5_Dt/test.out -J D5_Dt_${i} " 
     python D5_Dt_syntelog.py  D5_Dt/${i}  D5_Dt/${i}_out 
    "
done

for i in $(ls splitData); do
    bsub -q normal -n 1 -e splitData/test.err -o splitData/test.out -J splitData_${i} " 
     python syntelog_uniq.py  splitData/${i}  splitData/${i}_out 
    "
done

