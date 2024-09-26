
###
 # @Descripttion: 
 # @version: 
 # @Author: zpliu
 # @Date: 2023-10-08 10:51:26
 # @LastEditors: zpliu
 # @LastEditTime: 2024-09-26 09:59:33
 # @@param: 
### 

module load interproscan/5.25-64.0
#* interproscan
peptideFile=/cotton/public/public_data/GWAS_Fiber_ljy/HiFi_Gene_Annotation/A2_HiFi/J85/Final_Gene/J85.final.pep.fa
awk '{print $1}' ${peptideFile} |sed 's/\.1//g'|sed 's/\.$//g' > J85.new.pep.fa 

bsub -q normal -n 20 -e test.err -o test.out -J GO "
    interproscan.sh -appl Pfam-31.0 -i J85.new.pep.fa -iprlookup -goterms -cpu 20 -pa
"