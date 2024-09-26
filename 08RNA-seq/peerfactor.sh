###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2024-09-26 11:31:48
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 11:31:48
 # @@param: 
### 
for stage in 0DPA 4DPA 8DPA 12DPA 16DPA 20DPA; do
        bsub -q normal -n 1 -e test.err -o test.out -M 10G -J $stage "
        cat ../phenotype/${stage}_osca_efile.bed | cut -f1 --complement | awk 'NR>1{print \$0}' | sed 's/\t/,/g' >${stage}_peer_input.txt
        python2.7 peer_interface.py \
                -f  ${stage}_peer_input.txt \
                -n 20 -o peer_${stage}_factor_20
         "
done