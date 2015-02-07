_colors = { 'nw': 'r',
            'kT0.5': 'b',
            'kT0.4': 'c',
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

def color(method):
    if method in _colors:
        return _colors[method]
    return 'k' # everything else black

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
               'bubble_suppression': 'bubble suppression',
               'wang_landau': 'Wang-Landau',
               'vanilla_wang_landau': 'Vanilla Wang-Landau',
               'robustly_optimistic': 'robustly optimistic',
               'gaussian': 'gaussian method',
               'walker_optimization': 'walker optimization'}
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
    if method in ['wang_landau','vanilla_wang_landau','robustly_optimistic']:
        return color(method) + '+'
    if method in _colors:
        return color(method) + '.'
    return '.'
