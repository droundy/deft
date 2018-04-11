import numpy as np
import re
import sys

moves = re.compile(r"moves (\d+), energies_found (\d+), erange: (\d+) -> (\d+)")
print "# energies_found, t, ehi, elo"
for line in sys.stdin.readlines():
    match = moves.search(line)
    if match is not None:
        g = match.groups()
        print g[1], g[0], g[2], g[3]
    
    
    


