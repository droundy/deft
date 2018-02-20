import matplotlib.pyplot as plt

_colors = { 'sad': 'r',
            'sad-tm': 'r',
            'sad3': 'c',
            'sad3-tm': 'c',
            'sad3-s1': 'c',
            'sad3-s1-tm': 'c',
            'sad3-s2': 'c',
            'sad3-s2-tm': 'c',
            'wltmmc-0.8-1e-10': 'k',
            'wltmmc-0.8-1e-10-s1': 'k',
            'wltmmc': 'k',
            'vanilla_wang_landau': 'k',
            'tmmc': 'b',
            'tmi3': 'g',
            'toe3': 'tab:orange',
            'samc': 'm',
            'satmmc': 'y',
            'wltmmc-0.8-0.0001': 'tab:purple',
            'wltmmc-1-0.0001': 'tab:pink',
}

_linestyles = {
    'sad-tm': '--',
    'sad3-tm': '--',
    'sad3-s1': ':',
    'sad3-s1-tm': '--',
    'sad3-s2': '-.',
    'sad3-s2-tm': '--',
    'vanilla_wang_landau': '--',
}

_legend_order = [
    'sad', 'sad-tm',
    'sad3', 'sad3-tm',
    'sad3-s1',
    'sad3-s2',
    'tmmc',
    'wltmmc-0.8-1e-10',
    'wltmmc',
    'vanilla_wang_landau',
    'tmi3', 'toe3',
    'samc',
    'satmmc',
]

_legend_label = {
    'vanilla_wang_landau': 'WL',
    'sad': 'SAD',
    'sad-tm': 'SAD-TM',
    'wltmmc-0.8-1e-10': 'WLTMMC ($10^{-10}$ cutoff)',
    'tmmc': 'TMMC',
    'tmi3': 'TMI',
    'toe3': 'TOE',
    'samc': 'SAMC',
}

def fix_legend(method):
    if method in _legend_label:
        return _legend_label[method]
    return method

def legend_order(method):
    if method in _legend_order:
        return _legend_order.index(method)
    return len(_legend_order)

def color(m):
    if m in _colors:
        return _colors[m]
    return 'k'

def style_args(method):
    args = {'label': method}
    if method in _colors:
        args['color'] = _colors[method]
    if method in _linestyles:
        args['linestyle'] = _linestyles[method]
    return args

def plot(x, y, method=None):
    return plt.plot(x,y, **style_args(method))

def loglog(x, y, method=None):
    return plt.loglog(x,y, **style_args(method))

def legend():
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key = lambda t: legend_order(t[0])))
    labels = map(fix_legend, labels)
    plt.legend(handles, labels, loc='best')
