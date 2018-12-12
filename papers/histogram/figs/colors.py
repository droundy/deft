import matplotlib.pyplot as plt

_colors = { 'sad': 'r',
            'sad-tm': 'r',
            'sad3': 'c',
            'sad3-T13': 'c',
            'sad3-T13-slow': 'c',
            'sad-256': 'c',
            'sad3-test': 'c',
            'sad3-slow': 'c',
            'sad3-fast': 'c',
            'sad3-tm': 'c',
            'sad3-s1': 'c',
            'sad3-s1-tm': 'c',
            'sad3-s2': 'r',
            'sad3-s2-tm': 'c',
            'wltmmc-1e-10': 'k',
            'wltmmc-1e-10-slow': 'k',
            'wltmmc-1e-10-fast': 'k',
            'wltmmc-0.8-1e-10': 'k',
            'wltmmc-0.0001-T13': 'tab:purple',
            'wltmmc-0.0001-T13-slow': 'tab:purple',
            'wltmmc-1e-10-T13': 'k',
            'wltmmc-1e-10-T13-slow': 'k',
            'wltmmc': 'k',
            'vanilla_wang_landau': 'k',
            'vanilla_wang_landau-minE': 'k',
            'vanilla_wang_landau-slow': 'k',
            'vanilla_wang_landau-fast': 'k',
            'one_over_t_wang_landau-T13-t': 'xkcd:lime green',
            'vanilla_wang_landau-T13': 'k',
            'vanilla_wang_landau-T13-slow': 'k',
            'tmmc': 'b',
            'tmmc-slow': 'b',
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
            'samc-100000': 'gold',
            'samc-1e+06': 'g',
            'samc-1000-slow': 'r',
            'samc-10000-slow': 'tab:orange',
            'samc-100000-slow': 'gold',
            'samc-1000-fast': 'r',
            'samc-10000-fast': 'tab:orange',
            'samc-100000-fast': 'gold',
            'samc-1e+06': 'g',
            'samc-1e+07': 'midnightblue',
            'samc-1e+06-slow': 'g',
            'samc-1e+07-slow': 'midnightblue',
            'samc-1e+08': 'darkmagenta',
            'samc-1e+09': 'fuchsia',
            'samc-500-1e5': 'xkcd:red',
            'samc-1e5-256': 'gold',
            'samc-1e6-256': 'g',
            '1/sqrt(t)': 'xkcd:light gray'
}

_linestyles = {
    'sad-tm': '--',
    'sad3-tm': '--',
    'sad3-s1': ':',
    'sad3-s1-tm': '--',
    'sad3-s2': '-.',
    'sad3-s2-tm': '--',
    'sad3-slow': '-',
    'sad3-fast': '-.',
    'vanilla_wang_landau': '--',
    'vanilla_wang_landau-T13': '--',
    'vanilla_wang_landau-T13-slow': '--',
    'vanilla_wang_landau-minE': '--',
    'vanilla_wang_landau-slow': '--',
    'vanilla_wang_landau-fast': '--',
    'samc-1000': ':',
    'samc-10000': ':',
    'samc-100000': ':',
    'samc-1000-slow': ':',
    'samc-10000-slow': ':',
    'samc-100000-slow': ':',
    'samc-1e+06-slow': ':',
    'samc-1e+07-slow': ':',
    'samc-1000-fast': ':',
    'samc-10000-fast': ':',
    'samc-100000-fast': ':',
    'samc-1e+06': ':',
    'samc-1e+07': ':',
    'samc-1e+08': ':',
    '1/sqrt(t)': '-'
}

_legend_order = [
    'vanilla_wang_landau-s2',
    'sad-256',
    'sad3-T13',
    'sad3-T13-slow',
    'sad', 'sad-tm',
    'sad3-test', 'sad3-test',
    'sad3', 'sad3-tm',
    'sad3-slow', 'sad3-slow-tm',
    'sad3-fast', 'sad3-fast-tm',
    'sad3-s1',
    'sad3-s2',
    'vanilla_wang_landau',
    'vanilla_wang_landau-minE',
    'vanilla_wang_landau-slow',
    'vanilla_wang_landau-fast',
    'vanilla_wang_landau-T13',
    'vanilla_wang_landau-T13-slow',
    'tmmc',
    'tmmc-slow',
    'wltmmc-0.8-0.0001',
    'wltmmc-0.0001',
    'wltmmc-0.0001-slow',
    'wltmmc-0.0001-fast',
    'wltmmc-0.0001-T13',
    'wltmmc-0.0001-T13-slow',
    'wltmmc',
    'wltmmc-0.8-1e-10',
    'wltmmc-1e-10',
    'wltmmc-1e-10-slow',
    'wltmmc-1e-10-fast',
    'wltmmc-1e-10-T13',
    'wltmmc-1e-10-T13-slow',
    'tmi3', 'toe3',
    'samc',
    'samc-1000',
    'samc-10000',
    'samc-100000',
    'samc-1e+06',
    'samc-1e+07',
    'samc-1e+08',
    'samc-1e+09',
    'samc-1000-slow',
    'samc-10000-slow',
    'samc-100000-slow',
    'samc-1e+06-slow',
    'samc-1e+07-slow',
    'samc-1e+08-slow',
    'samc-1e+09-slow',
    'samc-1000-fast',
    'samc-10000-fast',
    'samc-100000-fast',
    'samc-500-1e5',
    'samc-1e5-256',
    'samc-1e6-256',
]

_legend_label = {
    'vanilla_wang_landau-T13': 'WL',
    'vanilla_wang_landau-T13-slow': 'WL',
    'vanilla_wang_landau': 'WL',
    'vanilla_wang_landau-minE': 'WL',
    'vanilla_wang_landau-slow': 'WL',
    'vanilla_wang_landau-fast': 'WL',
    'vanilla_wang_landau_C': 'WL-C',
    'vanilla_wang_landau-s2': 'WL-s2',
    'one_over_t_wang_landau-T13-t': '1/t-WL',
    'sad': 'SAD',
    'sad-256': 'SAD',
    'sad3-s2': 'SAD-s2',
    'sad-tm': 'SAD-TM',
    'sad3': 'SAD',
    'sad3-T13': 'SAD',
    'sad3-T13-slow': 'SAD',
    'sad3-slow': 'SAD',
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
    'wltmmc-1e-10-T13': 'WLTMMC ($10^{-10}$ cutoff)',
    'wltmmc-1e-10-T13-slow': 'WLTMMC ($10^{-10}$ cutoff)',
    'wltmmc-0.0001-T13': 'WLTMMC ($10^{-4}$ cutoff)',
    'wltmmc-0.0001-T13-slow': 'WLTMMC ($10^{-4}$ cutoff)',
    'tmmc': 'TMMC',
    'tmmc-slow': 'TMMC',
    'tmi3': 'TMI',
    'toe3': 'TOE',
    'samc': 'SAMC',
    'samc-1000': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1e+06': 'SAMC ($10^{6}$ $t_0$)',
    'samc-1e+07': 'SAMC ($10^{7}$ $t_0$)',
    'samc-1e+08': 'SAMC ($10^{8}$ $t_0$)',
    'samc-1e+09': 'SAMC ($10^{9}$ $t_0$)',
    'samc-1000-slow': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000-slow': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000-slow': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1e+06-slow': 'SAMC ($10^{6}$ $t_0$)',
    'samc-1e+07-slow': 'SAMC ($10^{7}$ $t_0$)',
    'samc-1e+08-slow': 'SAMC ($10^{8}$ $t_0$)',
    'samc-1e+09-slow': 'SAMC ($10^{9}$ $t_0$)',
    'samc-1000-fast': 'SAMC ($10^{3}$ $t_0$)',
    'samc-10000-fast': 'SAMC ($10^{4}$ $t_0$)',
    'samc-100000-fast': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1e5-256': 'SAMC ($10^{5}$ $t_0$)',
    'samc-1e6-256': 'SAMC ($10^{6}$ $t_0$)',
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
    args = {}
    if method != '1/sqrt(t)':
        args = {'label': method}
    if method in _colors:
        args['color'] = _colors[method]
    if method in _linestyles:
        args['linestyle'] = _linestyles[method]
    args['zorder'] = -legend_order(method)
    if method == '1/sqrt(t)':
        args['zorder'] = -50
    if method == '1/sqrt(t)':
        args['linewidth'] = 0.05
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
