# compute power of our statistic in windows, as compared to mixed background, neutral background, and selected background
import sys
import numpy as np
from collections import defaultdict

ll_data = defaultdict(dict)

# read in LL file, store data for each window
file = "/Users/telemacher/projects/2DSFS-scan/sim_likelihood_calcs/migRate."+sys.argv[1]+"."+sys.argv[2]+"_fst.txt"

fh = open(file, 'r')
for line in fh:
    break
for line in fh:
    data = line.strip().split()
    ll_data[float(data[3])][float(data[4])] = data[5]

# ignore windows used to compute BG distribution
# if window overlaps selected locus, put into selected fg or bg group
# if window overlaps neutral lcous only, put into neutral group

bg_NeutWeakSel = []
bg_LargeSel = [] 
fg_NeutWeakSel = []
fg_LargeSel = []

for window_start in ll_data:
    for window_end in ll_data[window_start]:
        if window_end < 500001:
            continue
        if window_start >= 500001 and window_end < 990000:
            fg_NeutWeakSel.append(ll_data[window_start][window_end])
        if window_start >= 990000 and window_end < 1000000:
            fg_LargeSel.append(ll_data[window_start][window_end])
        if window_start >= 1000001 and window_end < 1490000:
            bg_NeutWeakSel.append(ll_data[window_start][window_end])
        if window_start >= 1490000 and window_end < 1500000:
            bg_LargeSel.append(ll_data[window_start][window_end])

# compute percentiles of BG groups (0.95)
bg_NeutWeakSel.sort()
fg_NeutWeakSel.sort()
bg_LargeSel.sort()
fg_LargeSel.sort()

thresh_neut = bg_NeutWeakSel[int(0.995*len(bg_NeutWeakSel))]
thresh_sel = bg_LargeSel[int(0.995*len(bg_LargeSel))]

# compute power for each percentile group for each BG
power_neut = 0 
power_sel = 0 
for val in fg_LargeSel:
    if val > thresh_neut:
        power_neut += 1
    if val > thresh_sel:
        power_sel += 1

print(sys.argv[1],sys.argv[2],power_sel/len(fg_LargeSel), power_neut/len(fg_LargeSel))



