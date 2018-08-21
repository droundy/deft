import matplotlib.pyplot as plt

_colors = { 'sad': 'r',
            'sad-tm': 'r',
            'sad3': 'c',
            'sad3-slow': 'c',
            'sad3-fast': 'c',
            'sad3-tm': 'c',
            'sad3-s1': 'c',
            'sad3-s1-tm': 'c',
            'sad3-s2': 'c',
            'sad3-s2-tm': 'c',
            'wltmmc-1e-10': 'k',
            'wltmmc-1e-10-slow': 'k',
            'wltmmc-1e-10-fast': 'k',
            'wltmmc-0.8-1e-10': 'k',
            'wltmmc-0.8-1e-10-s1': 'k',
            'wltmmc': 'k',
            'vanilla_wang_landau': 'k',
            'vanilla_wang_landau-slow': 'k',
            'vanilla_wang_landau-fast': 'k',
            'tmmc': 'b',
            'tmi3': 'g',
            'toe3': 'tab:orange',
            'samc': 'm',
            'wltmmc-0.8-0.0001': 'tab:purple',
            'wltmmc-0.0001': 'tab:purple',
            'wltmmc-0.0001-slow': 'tab:purple',
            'wltmmc-0.0001-fast': 'tab:purple',
            'wltmmc-1-0.0001': 'tab:pink',
            'wltmmc': 'tab:purple',
            'samc-1000': 'r',
            'samc-10000': 'tab:orange',
            'samc-100000': 'y',
            'samc-1e+06': 'g',
            'samc-1000-slow': 'r',
            'samc-10000-slow': 'tab:orange',
            'samc-100000-slow': 'y',
            'samc-1000-fast': 'r',
            'samc-10000-fast': 'tab:orange',
            'samc-100000-fast': 'y',
            'samc-1e+06': 'g',
            'samc-1e+06': 'g',
            '1/sqrt(t)': 'r'
}

_linestyles = {
    'sad-tm': '--',
    'sad3-tm': '--',
    'sad3-s1': ':',
    'sad3-s1-tm': '--',
    'sad3-s2': '-.',
    'sad3-s2-tm': '--',
    'sad3-slow': '--',
    'sad3-fast': '-.',
    'vanilla_wang_landau': '--',
    'vanilla_wang_landau-slow': '--',
    'vanilla_wang_landau-fast': '--',
    'samc-1000': ':',
    'samc-10000': ':',
    'samc-100000': ':',
    'samc-1000-slow': ':',
    'samc-10000-slow': ':',
    'samc-100000-slow': ':',
    'samc-1000-fast': ':',
    'samc-10000-fast': ':',
    'samc-100000-fast': ':',
    'samc-1e+06': ':',
}

_legend_order = [
    'sad', 'sad-tm',
    'sad3', 'sad3-tm',
    'sad3-slow', 'sad3-slow-tm',
    'sad3-fast', 'sad3-fast-tm',
    'sad3-s1',
    'sad3-s2',
    'vanilla_wang_landau',
    'vanilla_wang_landau-slow',
    'vanilla_wang_landau-fast',
    'tmmc',
    'wltmmc-0.8-0.0001',
    'wltmmc-0.0001',
    'wltmmc-0.0001-slow',
    'wltmmc-0.0001-fast',
    'wltmmc',
    'wltmmc-0.8-1e-10',
    'wltmmc-1e-10',
    'wltmmc-1e-10-slow',
    'wltmmc-1e-10-fast',
    'tmi3', 'toe3',
    'samc',
    'samc-1000',
    'samc-10000',
    'samc-100000',
    'samc-1000-slow',
    'samc-10000-slow',
    'samc-100000-slow',
    'samc-1000-fast',
    'samc-10000-fast',
    'samc-100000-fast',
    'samc-1e+06',
]

_legend_label = {
    'vanilla_wang_landau': 'WL',
    'vanilla_wang_landau-slow': 'WL',
    'vanilla_wang_landau-fast': 'WL',
    'sad': 'SAD',
    'sad-tm': 'SAD-TM',
    'sad3': 'SAD',
    'sad3-slow': 'SAD-slow',
    'sad3-fast': 'SAD-fast',
    'wltmmc': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.8-0.0001': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.0001': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.0001-slow': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.0001-fast': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.8-1e-10': 'WLTMMC ($10^{-10}$ cutoff)',
    'wltmmc-1e-10': 'WLTMMC ($10^{-10}$ cutoff)',
    'wltmmc-1e-10-slow': 'WLTMMC ($10^{-10}$ cutoff)',
    'wltmmc-1e-10-fast': 'WLTMMC ($10^{-10}$ cutoff)',
    'tmmc': 'TMMC',
    'tmi3': 'TMI',
    'toe3': 'TOE',
    'samc': 'SAMC',
    'samc-1000': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1000-slow': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000-slow': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000-slow': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1000-fast': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000-fast': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000-fast': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1e+06': 'SAMC ($10^{6}$ $t_0$)',
    '1/sqrt(t)': r'$\frac{1}{\sqrt{t}}$'
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
    args['zorder'] = -legend_order(method)
    return args

def plot(x, y, method=None):
    return plt.plot(x,y, **style_args(method))

def loglog(x, y, method=None):
    return plt.loglog(x,y, **style_args(method))

def legend(loc='best'):
    handles, labels = plt.gca().get_legend_handles_labels()
    labels_and_handles = list(zip(*sorted(zip(labels, handles), key = lambda t: legend_order(t[0]))))
    if len(labels_and_handles) == 2:
      labels, handles = labels_and_handles
      labels = map(fix_legend, labels)
      plt.legend(handles, labels, loc=loc)
