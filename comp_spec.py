



# 
import sys
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
fh = open(sys.argv[1],'r')


back = defaultdict(dict)
for line in fh:
    break
for line in fh:
    data = line.strip().split()
    back[data[0]][data[1]] = [int(data[2]), float(data[3])]
fh.close() 


fh = open(sys.argv[2],'r')


tot = 0
fore = defaultdict(dict)
for line in fh:
    break
for line in fh:
    data = line.strip().split()
    fore[data[0]][data[1]] = [int(data[2]),float(data[3])]
    tot += float(data[2])

fh.close() 


for item in back:
    for item2 in back[item]:
        if back[item][item2][0] > 0:
            print(item,item2,(float(item2)+float(item))/(14+14), poisson.logpmf(k=fore[item][item2][0],mu=back[item][item2][1]*tot)*np.sign(fore[item][item2][1]-back[item][item2][1]), (fore[item][item2][1]-back[item][item2][1])/back[item][item2][0]**0.5)
        else:
            print(item,item2,(float(item2)+float(item))/(14+14), poisson.logpmf(k=fore[item][item2][0],mu=back[item][item2][1]*tot)*np.sign(fore[item][item2][1]-back[item][item2][1]),0)
            



