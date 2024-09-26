<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2024-09-26 10:24:48
 * @LastEditors: zpliu
 * @LastEditTime: 2024-09-26 12:12:47
 * @@param: 
-->
### `get_distances.py` 

> Calculate the maximum Jaccard similarity between cultivar group and the semi-wild group.


### Example
```bash
python get_distances.py -vcf 59AD1.MAF0.02.PAVs.clean.vcf -chr HC04_${ID} -species_file treeID.txt -queryspecies     Wild -targetspecies Cultivar -fai HC04_chr_adjust.fa.fai -w 1000000 >Wild_Cultivar.${ID}.similarity

```