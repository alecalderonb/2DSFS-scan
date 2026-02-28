import sys

fh = open(sys.argv[1],'r')

names = []

for line in fh:
    data = line.strip().split(',')
    for thing in data:
        thing = thing.split()
        if len(thing) < 2:
            thing.append('')
        names.append(thing)

names.sort(key=lambda x: x[1])
for name in names:
    print(name[0],name[1])

    

