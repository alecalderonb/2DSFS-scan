from twoDSFS import likelihoodCalcs
import bz2
import pickle
import sys

fn = sys.argv[1]
nSNPs = int(sys.argv[2])
inferencePipeline = likelihoodCalcs.LikelihoodInference_jointSFS(regen=False)

# read and write pickle data
inferencePipeline.make_data_dict_vcf(fn, popinfo_filename = 'data/popmap.txt',pickleWrite = True, outFile = 'data/pruned.'+str(nSNPs)+'.pkl.bz2')

# load ECB snp data 

with bz2.BZ2File('data/pruned.'+str(nSNPs)+'.pkl.bz2', 'rb') as file:
    ECB_wg_dict = pickle.load(file)

'''calculate 2D sfs'''
my2Dsfs = inferencePipeline.calculate_2d_sfs(ECB_wg_dict)

''' each chromosome as its own background '''
ECB_stats = inferencePipeline.T2D_scan(ECB_wg_dict, numSnps=nSNPs)

for wind in ECB_stats:
    print(wind,  ECB_stats[wind]['snp_count'],  ECB_stats[wind]['T2D'])
