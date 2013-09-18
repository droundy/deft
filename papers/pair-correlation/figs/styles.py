

color = { 'mc': 'k',
          'this-work': 'b',
          'fischer': 'g',
          'sokolowski': 'r'}

line = { 'mc': '.',
         'this-work': '-',
         'fischer': '-.',
         'sokolowski': '--'}

back_line = { 'mc': '.',
              'this-work': '-',
              'fischer': '-.',
              'sokolowski': '--'}

title = { 'mc': 'Monte Carlo',
          'this-work': 'this work',
          'fischer': 'Fischer',
          'sokolowski': 'Sokolowski'}

plot = {}
for k in color:
    plot[k] = color[k] + line[k]

plot_back = {}
for k in color:
    plot_back[k] = color[k] + back_line[k]
