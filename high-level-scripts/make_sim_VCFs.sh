#for migRate in {0.0,0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do python scripts/make_simulated_VCFs.py $migRate $gen; done; done
for migRate in {0.0,0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do python scripts/make_simulated_VCFs.py $migRate $gen True; done; done
