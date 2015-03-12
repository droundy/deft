_colors = { 'nw': 'r',
            'kT0.5': 'r',
            'kT0.4': 'c',
            'kT1': 'b',
            'kT2': 'c',
            'kT0.1': 'g',
            'tmmc': 'k',
            'oetmmc': 'b',
            'wang_landau': 'g',
            'vanilla_wang_landau': 'b',
            'simple_flat': 'r',
            'gaussian': 'm',
            'optimized_ensemble': 'y'}

def color(method):
    if method in _colors:
        return _colors[method]
    return 'b' # everything else blue

def line(method):
    lines = { 'nw': '--',
              'tmmc': ':'}
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
               'gaussian': 'gaussian method',
               'optimized_ensemble': 'optimized ensemble'}
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
    if method in ['wang_landau','vanilla_wang_landau','simple_flat']:
        return color(method) + '+'
    if method in _colors:
        return color(method) + '.'
    return '.'
