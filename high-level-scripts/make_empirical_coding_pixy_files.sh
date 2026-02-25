

#python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf.gz 100 > data/T2D.coding.100.txt 
#python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf.gz 250 > data/T2D.coding.250.txt 
#python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf.gz 500 > data/T2D.coding.500.txt 

# bgzip the data
#gunzip  data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf.gz
#cat data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf | awk '{if($1 !~ "##INFO=<ID=VDB,Number=1") print}' > data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.2.vcf 
#gzip data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.vcf
#bgzip -f data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.2.vcf	  

# use tabix to index the file
#tabix -f -p vcf data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.2.vcf.gz

# make the BED files
#for n in {100,250,500}; do cat data/T2D.coding.${n}.txt | awk '{print $1"\t"$2"\t"$3}' > data/T2D.coding.${n}.bed; done

#conda activate pixy

for n in {200,250}; do pixy --stats fst --vcf data/vcf/ECBAnnotated_biallelic.coding.pruned.r2_0.5.10kb.2.vcf.gz --populations  data/popmap.noUnder.txt --bed_file data/T2D.coding.${n}.bed --output_folder data/ECB_pixy  --output_prefix ECB.coding.${n} --bypass_invariant_check; done

#conda deactivate

