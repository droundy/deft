def density_color(rd):
    dc = { '0.1': 'k',
           '0.2': 'y',
           '0.3': 'g',
           '0.4': 'r',
           '0.5': 'm',
           '0.6': 'c',
           '0.7': 'b',
           '0.8': 'k',
           '0.9': 'y',
           '1.0': 'r',
        }
    key = '%.5g' % rd
    if key in dc:
        return dc[key]
    return 'k'

color = { 100.0: 'xkcd:dark blue',
          20.0: 'xkcd:forest green',
          10.0: 'c',
          5.0: 'xkcd:orange',
          3.0: 'xkcd:lavender',
          2.5: 'xkcd:hot pink',
          2.0: 'g',
          1.5: 'xkcd:magenta',
          1.0: 'b',
          0.6: 'm',
          0.5: 'r',
          0.4: 'y',
          0.2: 'k',
          0.1: 'r',
          0.03: 'y',
          0.01: 'c',
          0.001: 'b',
          0.0001: 'c',
          0.00001: 'g',
          0.0: 'k' }
line = { 'wcadft': '-',
         'wcamc' :'--'}

def new_dft_linestyle():
    return '-'

def mcwca_linestyle():
    return '--'

def other_mcwca_linestyle():
    return '.'

def bh_dft_linestyle():
    return ':'

def color_from_kT(kT):
    if kT in color:
        return color[kT]
    return 'xkcd:vomit'