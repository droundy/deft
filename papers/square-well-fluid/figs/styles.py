color = { 'mc': 'k',
          'this-work': 'b',
          'this-work-mc': 'c',
          'this-work-short': 'b',
          'fischer': 'g',
          'sokolowski': 'r'}

line = { 'mc': '.',
         'this-work': '-',
         'this-work-mc': ':',
         'this-work-short': '-',
         'fischer': '-.',
         'sokolowski': '--'}

forward_line = { 'mc': '.',
                 'this-work': '>',
                 'this-work-mc': '>',
                 'this-work-short': '>',
                 'fischer': '-.',
                 'sokolowski': '>'}

back_line = { 'mc': '.',
              'this-work': '<',
              'this-work-mc': '<',
              'this-work-short': '<',
              'fischer': '-.',
              'sokolowski': '<'}

title = { 'mc': 'Monte Carlo',
          'this-work': 'CVA-S',
          'this-work-mc': 'CVA-old',
          'this-work-short': 'CVA',
          'fischer': 'Fischer',
          'sokolowski': 'Sokolowski'}

start = { 'this-work-short' : 0,
          'this-work-mc' : 1/6.,
          'sokolowski' : 1/3.}

plot = {}
for k in color:
    plot[k] = color[k] + line[k]

plot_back = {}
for k in color:
    plot_back[k] = color[k] + back_line[k]

plot_forward = {}
for k in color:
    plot_forward[k] = color[k] + forward_line[k]

plots = ['mc', 'sokolowski', 'fischer', 'this-work-short']
oriented_plots = ['sokolowski', 'this-work-short']

short_range = 4.0
