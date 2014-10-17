
color = { 100.0: 'k',
          10.0: 'c',
          2.0: 'g',
          1.0: 'k',
          0.6: 'r',
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
new_dft_code = {}
dftwca = {}
mc = {}
mcwca = {}

new_dft_code = {}
for k in color:
    dftwca[k] = color[k] + ':'
    new_dft_code[k] = color[k] + ':.'
    dft[k] = color[k] + '--'
    coarsedft[k] = color[k] + '-'
    mc[k] = color[k] + '.'
    mcwca[k] = color[k] + '-'
