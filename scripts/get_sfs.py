#

import sys
import gzip as gz
from collections import defaultdict

def get_pop1_freq(genos):

    freq = 0
    for geno in genos:
        alleles = geno.split('|')
        alleles = [int(alleles[0]),int(alleles[1])]
        freq += alleles[0]
        freq += alleles[1]

    return freq
 

# open vcf and read in data for snp_type of interest

sfs = defaultdict(dict)
annot = {}
annot[(990000,1000000)] = "fg"
annot[(1490000,1500000)] = "bg"

for gen in (3000,3180,3300):
    for coords in [(990000,1000000),(1490000,1500000)]:
        sfs[annot[coords]][gen] = [0 for i in range(20)]  
        with gz.open("slim_sims/LHU_sims/simulations_LHU/vcfs_dadiparams/migRate.0.01."+str(gen)+".vcf.gz", 'rt') as f:
            for line in f:
                if line[0] == "#":
                    continue
                data = line.strip().split()
                annots = data[7].split(';')
                mybreak = False
                for myannot in annots:
                    keyVal = myannot.split('=')
                    if keyVal[0] == "S":
                        if float(keyVal[1]) == 0.:
                            mybreak = True
                if mybreak:
                    continue            
                if int(data[1]) >= coords[0] and int(data[1]) <= coords[1]:
                    freq = get_pop1_freq(data[9:19]) 
                    if freq < 20 and freq > 0:
                        sfs[annot[coords]][gen][freq] +=1

for annot in sfs:
    for gen in sfs[annot]:
        tot = sum(sfs[annot][gen])
        for i in range(len(sfs[annot][gen])):
            if i > 0:
                print(annot, gen, i, sfs[annot][gen][i]/tot)
