color = { 'nw': 'r',
          'kT1': 'b',
          'kT2': 'c',
          'kT0.1': 'g',
          'tmmc': 'k',
          'bubble_suppression': 'k',
          'wang_landau': 'g',
          'vanilla_wang_landau': 'b',
          'robustly_optimistic': 'r',
          'gaussian': 'm',
          'walker_optimization': 'y'}

line = { 'nw': '--',
         'kT1': '--',
         'kT2': ':',
         'kT0.1': '-.',
         'tmmc': '-',
         'bubble_suppression': '-',
         'wang_landau': '-',
         'vanilla_wang_landau': '-',
         'robustly_optimistic': '-',
         'gaussian': '-',
         'walker_optimization': '-'}

title = { 'nw': '$kT/\epsilon = \infty$ sim.',
          'kT1': '$kT/\epsilon = 1$ sim.',
          'kT2': '$kT/\epsilon = 2$ sim.',
          'kT0.1': '$kT/\epsilon = 0.1$ sim.',
          'tmmc': 'tmmc',
          'bubble_suppression': 'bubble suppression',
          'wang_landau': 'Wang-Landau',
          'vanilla_wang_landau': 'Vanilla Wang-Landau',
          'robustly_optimistic': 'robustly optimistic',
          'gaussian': 'gaussian method',
          'walker_optimization': 'walker optimization'}

plot = {}
for k in color:
    plot[k] = color[k] + line[k]

dots = {}
for k in color:
    dots[k] = color[k] + '.'
    if k in ['wang_landau','vanilla_wang_landau','robustly_optimistic']:
        dots[k] = color[k] + '+'
