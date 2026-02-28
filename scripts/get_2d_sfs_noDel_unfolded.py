import sys
import gzip as gz
from collections import defaultdict

output_file = sys.argv[1]

def get_pop_freq(genos):
    freq = 0
    for geno in genos:
        a1, a2 = geno.split('|')
        freq += int(a1)
        freq += int(a2)
    return freq

sfs = defaultdict(dict)
annot = {}
annot[(990000,1000000)] = "fg"
annot[(1490000,1500000)] = "bg"
#annot[(1000000,1490000)] = "bg" # null locus with low-effect alleles

for gen in (3000, 3180, 3300, 3600):
    for coords in [(990000,1000000),(1490000,1500000)]:
    #for coords in [(990000,1000000),(1000000,1490000)]:
          
        # 2D SFS: pop1 x pop2
        sfs[annot[coords]][gen] = [[0 for j in range(21)] for i in range(21)]
        
        with gz.open(f"/cluster/tufts/uricchiolab/ECB/simulations_LHU/scripts_AC/migRate.0.01.{gen}.noDel.vcf.gz", "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                
                data = line.strip().split()
                annots = data[7].split(';')
                skip = False
                for a in annots:
                    keyVal = a.split('=')
                    if keyVal[0] == "S": 
                        if float(keyVal[1]) == 0.:
                            skip = True
                if skip:
                    continue
                
                pos = int(data[1])
                if coords[0] <= pos <= coords[1]:
                    
                    freq_p1 = get_pop_freq(data[9:19])
                    freq_p2 = get_pop_freq(data[19:29])
                    
                    sfs[annot[coords]][gen][freq_p1][freq_p2] += 1
                        
                        
with open(output_file, "w") as out:
    out.write("region\tgeneration\tfreq_p1\tfreq_p2\tcount\tdensity\n")

    for region in sfs:
        for gen in sfs[region]:

            total = sum(
                sfs[region][gen][i][j]
                for i in range(21)
                for j in range(21)
            )

            for i in range(21):
                for j in range(21):
                    if i > 0 and j > 0: #and sfs[region][gen][i][j] > 0:
                        count = sfs[region][gen][i][j]
                        density = sfs[region][gen][i][j] / total
                        out.write(
                            f"{region}\t{gen}\t{i}\t{j}\t{count}\t{density}\n"
                        )    
                