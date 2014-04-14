color = { '-nw': 'r',
          '-kT1': 'b',
          '-kT2': 'c',
          '-kT0.1': 'g',
          '-flat': 'k'}

line = { '-nw': '--',
         '-kT1': '--',
         '-kT2': ':',
         '-kT0.1': '-.',
         '-flat': '-'}

title = { '-nw': '$kT/\epsilon = \infty$ sim.',
          '-kT1': '$kT/\epsilon = 1$ sim.',
          '-kT2': '$kT/\epsilon = 2$ sim.',
          '-kT0.1': '$kT/\epsilon = 0.1$ sim.',
          '-flat': 'flat histogram'}

plot = {}
for k in color:
    plot[k] = color[k] + line[k]

dots = {}
for k in color:
    dots[k] = color[k] + '.'
