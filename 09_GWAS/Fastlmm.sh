module load Fastlmm/v0.2.32
genotype='Ga_pan'
covarFile='K2_215Ga_StrucMatrix.txt'

phenotype='All_trait_blup_normal.txt'
outPath='./'

less ${phenotype} |
    head -1 | sed 's/\t/\n/g' |
    awk 'NR>=3{print $0}' |
    awk '{print NR"\t"$0}' |
    while read traitId trait; do
        fastlmmc -bfile ${genotype} -bfilesim ${genotype} \
            -pheno ${phenotype} -covar ${covarFile} \
            -mpheno ${traitId} \
            -maxThreads 1 \
            -out ${outPath}/blup_SV_${trait}.txt
    done
	