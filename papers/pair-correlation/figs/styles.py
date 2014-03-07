

color = { 'mc': 'k',
          'this-work': 'b',
          'this-work-mc': 'k',
          'fischer': 'g',
          'sokolowski': 'r'}

line = { 'mc': '.',
         'this-work': '-',
         'this-work-mc': ':',
         'fischer': '-.',
         'sokolowski': '--'}

forward_line = { 'mc': '.',
                 'this-work': '->',
                 'this-work-mc': ':>',
                 'fischer': '-.',
                 'sokolowski': '-->'}

back_line = { 'mc': '.',
              'this-work': '-<',
              'this-work-mc': ':<',
              'fischer': '-.',
              'sokolowski': '--<'}

title = { 'mc': 'Monte Carlo',
          'this-work': 'this work',
          'this-work-mc': 'this work MC',
          'fischer': 'Fischer',
          'sokolowski': 'Sokolowski'}

plot = {}
for k in color:
    plot[k] = color[k] + line[k]

plot_back = {}
for k in color:
    plot_back[k] = color[k] + back_line[k]

plot_forward = {}
for k in color:
    plot_forward[k] = color[k] + forward_line[k]
