hist_methods = ['wang_landau','vanilla_wang_landau','simple_flat','tmmc','oetmmc']

_colors = { 'nw': 'r',
            'kT0.4': 'c',
            'kT0.5': 'y',
            'kT1': 'g',
            'wang_landau': 'g',
            'simple_flat': 'r',
            'tmmc': 'k',
            'oetmmc': 'm',
            'wang_landau_oe': 'g',
            'simple_flat_oe': 'r',
            'tmmc_oe': 'k',
            'oetmmc_oe': 'm',
            'vanilla_wang_landau': 'y',
            'tmmc-golden': 'b',
            'cfw': 'c'}

def color(method):
    if method in _colors:
        return _colors[method]
    return 'b' # everything else blue

def line(method):
    if method in ['nw'] + [ h+'_oe' for h in hist_methods ] or method[:2] == 'kT':
        return '--'
    return '-'

def title(method):
    titles = { 'nw': '$kT/\epsilon = \infty$ sim.',
               'wang_landau': 'Wang-Landau',
               'simple_flat': 'Simple Flat',
               'tmmc': 'TMMC',
               'oetmmc': 'OETMMC',
               'wang_landau_oe': 'Wang-Landau, OE',
               'simple_flat_oe': 'Simple Flat, OE',
               'tmmc_oe': 'TMMC, OE',
               'oetmmc_oe': 'OETMMC, OE',
               'vanilla_wang_landau': 'Vanilla Wang-Landau',
               'tmmc-golden': 'TMMC Golden',
               'cfw': 'Converged Flat'}
    if method in titles:
        return titles[method]
    if method[:2] == 'kT':
        return '$kT/\epsilon = %g$ sim.' % float(method[2:])
    return 'unrecognized method'

def plot(method):
    if method in _colors:
        return color(method) + line(method)
    return line(method)

def dots(method):
    if method in hist_methods:
        return color(method) + '+'
    if method in [ h+'_oe' for h in hist_methods ]:
        return color(method) + 'x'
    if method in _colors:
        return color(method) + '.'
    return '.'
