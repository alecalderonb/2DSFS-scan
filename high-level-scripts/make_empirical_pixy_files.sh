# bgzip the data
#bgzip -f data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf	  

# use tabix to index the file
#tabix -f -p vcf data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz

# make the BED files
#for n in {250,500,1000,2000}; do cat data/T2D.${n}.txt | awk '{print $1"\t"$2"\t"$3}' > data/T2D.${n}.bed; done

#for n in {250,500,1000,2000}; do pixy --stats fst --vcf data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz --populations  data/popmap.noUnder.txt --bed_file data/T2D.${n}.bed --output_folder data/ECB_pixy  --output_prefix ECB.${n} --bypass_invariant_check; done


