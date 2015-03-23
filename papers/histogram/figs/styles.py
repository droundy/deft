_colors = { 'nw': 'r',
            'kT0.4': 'c',
            'kT0.5': 'y',
            'kT1': 'g',
            'tmmc': 'k',
            'oetmmc': 'm',
            'wang_landau': 'g',
            'vanilla_wang_landau': 'y',
            'simple_flat': 'r',
            'tmmc-golden': 'b'}

def color(method):
    if method in _colors:
        return _colors[method]
    return 'b' # everything else blue

def line(method):
    lines = { 'nw': '--'}
    if method[:2] == 'kT':
        return '--'
    if method in lines:
        return lines[method]
    return '-'

def title(method):
    titles = { 'nw': '$kT/\epsilon = \infty$ sim.',
               'tmmc': 'tmmc',
               'oetmmc': 'oetmmc',
               'wang_landau': 'Wang-Landau',
               'vanilla_wang_landau': 'Vanilla Wang-Landau',
               'simple_flat': 'simple flat',
               'tmmc-golden': 'tmmc golden'}
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
    if method in ['wang_landau','vanilla_wang_landau','simple_flat','tmmc','oetmmc']:
        return color(method) + '+'
    if method in _colors:
        return color(method) + '.'
    return '.'
