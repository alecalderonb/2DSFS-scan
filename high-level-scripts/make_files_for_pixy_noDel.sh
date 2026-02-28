# Make sorted vcf and zip with bgzip 

migRate=$1
gen=$2

gunzip -c slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams_noDel/migRate.${migRate}.${gen}.vcf.gz  | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > slim_vcfs/migRate.noDel.${migRate}.${gen}.vcf
bgzip  -f slim_vcfs/migRate.noDel.${migRate}.${gen}.vcf

# run our likelihood calculatiopn on this sorted vcf
python scripts/run_analyses.py False slim_vcfs/migRate.noDel.${migRate}.${gen}.vcf.gz > sim_likelihood_calcs/migRate.noDel.${migRate}.${gen}.temp.ll

# get rid of annoying warning message
cat sim_likelihood_calcs/migRate.noDel.${migRate}.${gen}.temp.ll | awk '{if(NR> 2) print}' > sim_likelihood_calcs/migRate.noDel.${migRate}.${gen}.ll

# index sorted vcf with tabix
tabix -f -p vcf slim_vcfs/migRate.noDel.${migRate}.${gen}.vcf.gz

# make bed file with appropriate format
cat sim_likelihood_calcs/migRate.noDel.${migRate}.${gen}.ll | awk '{if(int($2) == int($3)) {print $1"\t"int($2)"\t"int($3)+1} else {print $1"\t"int($2)"\t"int($3)}}' > slim_vcfs/migRate.noDel.${migRate}.${gen}.bed

# run pixy 
pixy --stats fst --vcf slim_vcfs/migRate.noDel.${migRate}.${gen}.vcf.gz --populations popmap.pixy.txt --bed_file slim_vcfs/migRate.noDel.${migRate}.${gen}.bed --output_folder sim_likelihood_calcs  --output_prefix migRate.noDel.${migRate}.${gen} --bypass_invariant_check



