<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2024-09-26 10:29:57
 * @LastEditors: zpliu
 * @LastEditTime: 2024-09-26 11:26:21
 * @@param: 
-->

### PAV hotspots

+ `PAV_hotspots/getData.ipynb` Get PAV hotspots based on GenomicRanges package in R 



### SYN and HYD

+ `SYN_HYD/A2_At_AllAnnotate.txt.gz` SYN and HYD between A2 and At 
+ `SYN_HYD/getData.ipynb` Annotated region for SYN and HYD
+ `SYN_HYD/A2_At_SYN_length` script for figure
+ `SYN_HYD/Circos` script for figure
+ `plot.ipynb` script for figure

```bash
#TODO get PSI region
module load cactus/2.5.1 
        halSynteny {input.halFile}  --maxAnchorDistance 5000 \
        --minBlockSize 5000 \
        --queryGenome HC04_At \
        --targetGenome DC085 \
        syn01_HC04_At_query_J85_target.psl

#TODO get SYN region in DC085 and HC04_At 
awk 'NR>1{print $14"\t"$16"\t"$17"\t"$10"\t"$12"\t"$13"\t"$9}' syn01_HC04_At_query_J85_target.psl  >J85_vs_HC04_synteny.txt

#TODO get HYD region in DC085 and HC04_At 
python SYN_HYD/intersectedBed_in_Python/main.py
```

### PAV diversity

+ `PAV_diversity/A2_PAV_pi.py` 
+ `PAV_diversity/At_PAV_pi.py` 
+ `PAV_diversity/Dt_PAV_pi.py` 


### Syntelog Gene analysis
