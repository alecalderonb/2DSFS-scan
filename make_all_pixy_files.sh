#for migRate in {0.0,0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do ./make_files_for_pixy.sh $migRate $gen; done; done
#for migRate in {0.01,0.05,0.1,0.2}; do for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do ./make_files_for_pixy_noDel.sh $migRate $gen; done; done
for gen in {3000,3060,3120,3180,3240,3300,3360,3420,3480,3540,3600}; do ./make_files_for_pixy_noDel.sh 0.2 $gen; done
for gen in {3420,3480,3540,3600}; do for migRate in {0.0,0.1}; do ./make_files_for_pixy_noDel.sh $migRate $gen; done; done
