python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz 250 > data/T2D.250.txt
python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz 500 > data/T2D.500.txt 
python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz 1000 > data/T2D.1000.txt
python scripts/run_empirical_analyses.py data/vcf/ECBAnnotated_biallelic.pruned.r2_0.5.10kb.vcf.gz 2000 > data/T2D.2000.txt

python combine_data.py data/ECB_pixy/ECB.250_fst.txt data/T2D.250.txt > data/windows.250.txt
python combine_data.py data/ECB_pixy/ECB.500_fst.txt data/T2D.500.txt > data/windows.500.txt
python combine_data.py data/ECB_pixy/ECB.1000_fst.txt data/T2D.1000.txt > data/windows.1000.txt
python combine_data.py data/ECB_pixy/ECB.2000_fst.txt data/T2D.2000.txt > data/windows.2000.txt
