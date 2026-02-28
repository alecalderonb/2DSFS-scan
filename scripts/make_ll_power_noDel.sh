rm  ../data/power.noDel.ll
rm   ../data/power.noDel.fst

for migRate in {0.0,0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do python  computePowerLL.noDel.py $migRate $gen >> ../data/power.noDel.ll ; done; done
for migRate in {0.0,0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do python  computePowerFst.noDel.py $migRate $gen >> ../data/power.noDel.fst ; done; done
