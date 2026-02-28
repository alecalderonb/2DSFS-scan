from twoDSFS import likelihoodCalcs
import bz2
import pickle
import sys

migRate = sys.argv[1]
gen = sys.argv[2]
noDel = False 
if len(sys.argv) > 3:
    noDel = bool(sys.argv[3])

inferencePipeline = likelihoodCalcs.LikelihoodInference_jointSFS()

# read and write pickle data
inferencePipeline.make_composite_VCF(migRate, gen, noDel = noDel)

"""

make_data_dict_vcf("slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams/migRate.0.03/iter_1/gen.8000.1.vcf.gz", "data/popmapSims.txt", pickleWrite = True, outFile = 'data/sim.pkl.bz2')

# load ECB snp data 

with bz2.BZ2File('data/sim.pkl.bz2', 'rb') as file:
    ECB_wg_dict = pickle.load(file)

'''calculate 2D sfs'''
my2Dsfs = inferencePipeline.calculate_2d_sfs(ECB_wg_dict)

''' each chromosome as its own background '''
ECB_stats_100kb = inferencePipeline.T2D_scan(ECB_wg_dict, my2Dsfs, numSnps=100)

for wind in ECB_stats_100kb:
    print(wind,  ECB_stats_100kb[wind]['snp_count'],  ECB_stats_100kb[wind]['T2D'])
"""

"""
plot_manhattan(ECB_stats_500kb, chr_ids, 'T1D_pop1', "univoltine T1D - 500kb windows - indep background")
plot_manhattan(ECB_stats_500kb, chr_ids, 'T1D_pop2', "bivoltine T1D - 500kb windows - indep background")
plot_manhattan(ECB_stats_500kb, chr_ids, 'T2D', title=None, threshold=5, ylim=(0,10000))
plot_manhattan(ECB_stats_500kb, chr_ids, 'new_term_pop1', "univoltine new_term - 500kb windows - indep background", ylim=(0,10000))
plot_manhattan(ECB_stats_500kb, chr_ids, 'new_term_pop2', "bivoltine new_term - 500kb windows - indep background", ylim=(0,10000))
plot_manhattan(ECB_stats_500kb, chr_ids, 'T2D_diff', "T2D - (T1Dpop1 + T1Dpop2)/2")

ECB_stats_20kb = inferencePipeline.combined_scan(ECB_wg_dict, 20000)

plot_manhattan(ECB_stats_20kb, chr_ids, 'T1D_pop1', "univoltine T1D - 20kb windows - indep background", ylim=(0,2200))
plot_manhattan(ECB_stats_20kb, chr_ids, 'T1D_pop2', "bivoltine T1D - 20kb windows - indep background", ylim=(0,2200))
plot_manhattan(ECB_stats_20kb, chr_ids, 'T2D', "T2D - 20kb windows - indep background")
plot_manhattan(ECB_stats_20kb, chr_ids, 'new_term_pop1', "univoltine new_term - 20kb windows - indep background", ylim=(0,2200))
plot_manhattan(ECB_stats_20kb, chr_ids, 'new_term_pop2', "bivoltine new_term - 20kb windows - indep background", ylim=(0,2200))
plot_manhattan(ECB_stats_20kb, chr_ids, 'T2D_diff', "T2D - (T1Dpop1 + T1Dpop2)/2 - 20kb windows")

save_csv_stats(ECB_stats_500kb, "/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/data/ECBstats_500kb.csv")

save_csv_stats(ECB_stats_20kb, "/Users/marlonalejandrocalderonbalcazar/Desktop/ECB/2DSFS_scan/data/ECBstats_20kb.csv")

ECB_stats_100kb = inferencePipeline.combined_scan(ECB_wg_dict, 100000)
plot_manhattan(ECB_stats_100kb, chr_ids, 'T1D_pop1', "univoltine T1D - 100kb windows")
plot_manhattan(ECB_stats_100kb, chr_ids, 'T1D_pop2', "bivoltine T1D - 100kb windows")
plot_manhattan(ECB_stats_100kb, chr_ids, 'T2D', "T2D - 100kb windows")
plot_manhattan(ECB_stats_100kb, chr_ids, 'new_term_pop1', "univoltine new_term - 100kb windows")
plot_manhattan(ECB_stats_100kb, chr_ids, 'new_term_pop2', "bivoltine new_term - 100kb windows")

''' specifying which chromosome to use as background '''
ECB_bgchr1_500kb = inferencePipeline.scan_chooseChr(ECB_wg_dict, 500000, 'NC_087088.1')
plot_manhattan(ECB_bgchr1_500kb, chr_ids, 'T1D_pop1', "univoltine T1D - chr1 background - 500kb windows")
plot_manhattan(ECB_bgchr1_500kb, chr_ids, 'T1D_pop2', "bivoltine T1D - chr1 background - 500kb windows")
plot_manhattan(ECB_bgchr1_500kb, chr_ids, 'T2D', "T2D - chr1 background - 500kb windows", threshold=5, ylim=(0,26000))
plot_manhattan(ECB_bgchr1_500kb, chr_ids, 'new_term_pop1', "univoltine new_term - chr1 background - 500kb windows")
"""
