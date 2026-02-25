import sys
from collections import defaultdict

fh = open(sys.argv[1], 'r')

fst = defaultdict(dict)

for line in fh:
    break

for line in fh:
    data = line.strip().split()
    fst[data[3]][data[4]] = [data[5],data[2]]
fh.close()

fh = open(sys.argv[2], 'r')

T2D = defaultdict(dict)

for line in fh:
    data = line.strip().split()
    T2D[data[1]][data[2]] = [data[4], data[3]]

for item in T2D:
    if item in fst:
        for item2 in T2D[item]:
            if item2 in fst[item]:
                print ( fst[item][item2][1], item, item2,  T2D[item][item2][0], fst[item][item2][0], T2D[item][item2][1])






    


