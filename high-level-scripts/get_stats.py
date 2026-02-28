
import sys

fh = open(sys.argv[1], 'r')
data = []

for line in fh:
    data.append([float(i) for i in line.strip().split()[2:4]])

fh.close()

fh = open(sys.argv[2], 'r')
for line in fh:
    break

T2Ds = {}
FSTs = {}

for line in fh:
    myline = line.strip().split()
    T2D = float(myline[4])
    FST = float(myline[5])
    T2Ds[T2D] = 0.
    FSTs[FST] = 0.
    

for val in data:
    i = 0
    for T2D in T2Ds:
        if T2D < val[0]:     
            i += 1
    pvalT2D = 1-i/(len(T2Ds)+1)
    i = 0
    for FST in FSTs:
        if FST < val[1]:     
            i += 1
    pvalFST = 1-i/(len(T2Ds)+1)

    print(val[0], pvalT2D, val[1], pvalFST)

