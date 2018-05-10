import numpy as np
import re
import sys

filename = str(sys.argv[1])

move = re.compile(r"We reached WL flatness (.*)")
wl = []
moves = []
for line in sys.stdin.readlines():
    match = move.search(line)
    if match is not None:
        g = match.groups()
        G = ''.join(map(str, g))
        movestemp = G.split(' ')[0]
        moves.append(int(movestemp.split('(')[-1]))
        wltemp = G.split('r')[-1]
        wl.append(np.double(wltemp.split(')')[0]))

np.savetxt('data/gamma/%s/wl.txt'%(filename), np.c_[moves,wl], 
delimiter = '\t', header = 'moves\t wl_factor')




