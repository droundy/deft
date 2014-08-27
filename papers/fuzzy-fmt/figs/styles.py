
color = { 0.6: 'r',
          0.4: 'y',
          0.2: 'k',
          0.1: 'r',
          0.03: 'y',
          0.01: 'm',
          0.001: 'b',
          0.0001: 'c',
          0.00001: 'g',
          0.0: 'k' }

coarsedft = {}
dft = {}
newdft = {}
mc = {}
mcljr = {}
for k in color:
    newdft[k] = color[k] + ':'
    dft[k] = color[k] + '--'
    coarsedft[k] = color[k] + '-'
    mc[k] = color[k] + '.'
    mcljr[k] = color[k] + '+'
