
color = { 0.0: 'k',
          0.1: 'r',
          0.03: 'y',
          0.01: 'm',
          0.001: 'b',
          0.0001: 'c',
          0.00001: 'g' }

coarsedft = {}
dft = {}
newdft = {}
mc = {}
for k in color:
    newdft[k] = color[k] + ':'
    dft[k] = color[k] + '--'
    coarsedft[k] = color[k] + '.--'
    mc[k] = color[k] + '-'
