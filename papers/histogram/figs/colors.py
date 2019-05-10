import matplotlib.pyplot as plt

_colors = { 'SAD': 'tab:orange',
            'WL': 'tab:green',
            '$1/t$-WL': 'tab:blue',
            'SAMC ($t_0 =10^{3}$)': 'r',
            'SAMC ($t_0 =10^{4}$)': 'xkcd:orange',
            'SAMC ($t_0 =10^{5}$)': 'xkcd:gold',
            'SAMC ($t_0 =10^{6}$)': 'g',
            'SAMC ($t_0 =10^{7}$)': 'midnightblue',
            'SAMC ($t_0 =10^{8}$)': 'darkmagenta',
            '1/sqrt(t)': 'xkcd:dark gray',
            r'$\frac{1}{\sqrt{t}}$': '#eeeeee',
}

_linestyles = {
    'SAMC ($t_0 =10^{3}$)': ':',
    'SAMC ($t_0 =10^{4}$)': ':',
    'SAMC ($t_0 =10^{5}$)': ':',
    'SAMC ($t_0 =10^{6}$)': ':',
    'SAMC ($t_0 =10^{7}$)': ':',
    'SAMC ($t_0 =10^{8}$)': ':',
    'WL': '--',
    '$1/t$-WL': '-.',
    r'$\frac{1}{\sqrt{t}}$': '-'
}

_legend_order = [
    'sad-256',
    'sad-50',
    'sad-50-slow',
    'wl-256',
    'wl-50',
    'wl-50-slow',
    'wl-inv-t-256',
    'wl-inv-t-50',
    'wl-inv-t-50-slow',
    'samc-1e3-256',
    'samc-1e4-256',
    'samc-1e5-256',
    'samc-1e6-256',
    'samc-1e7-256',
    'samc-1e8-256',
    'samc-1e3-50',
    'samc-1e4-50',
    'samc-1e5-50',
    'samc-1e6-50',
    'samc-1e7-50',
    'samc-1e3-50-slow',
    'samc-1e4-50-slow',
    'samc-1e5-50-slow',
    'samc-1e6-50-slow',
    'samc-1e7-50-slow',
    '1/sqrt(t)',
]

_legend_label = {
    'sad-256': 'SAD',
    'sad-50': 'SAD',
    'sad-50-slow': 'SAD',
    'wl-256': 'WL',
    'wl-50': 'WL',
    'wl-50-slow': 'WL',
    'wl-inv-t-256': '$1/t$-WL',
    'wl-inv-t-50': '$1/t$-WL',
    'wl-inv-t-50-slow': '$1/t$-WL',
    'samc-1e3-256': 'SAMC ($t_0 =10^{3}$)',
    'samc-1e4-256': 'SAMC ($t_0 =10^{4}$)',
    'samc-1e5-256': 'SAMC ($t_0 =10^{5}$)',
    'samc-1e6-256': 'SAMC ($t_0 =10^{6}$)',
    'samc-1e7-256': 'SAMC ($t_0 =10^{7}$)',
    'samc-1e8-256': 'SAMC ($t_0 =10^{8}$)',
    'samc-1e3-50': 'SAMC ($t_0 =10^{3}$)',
    'samc-1e4-50': 'SAMC ($t_0 =10^{4}$)',
    'samc-1e5-50': 'SAMC ($t_0 =10^{5}$)',
    'samc-1e6-50': 'SAMC ($t_0 =10^{6}$)',
    'samc-1e7-50': 'SAMC ($t_0 =10^{7}$)',
    'samc-1e3-50-slow': 'SAMC ($t_0 =10^{3}$)',
    'samc-1e4-50-slow': 'SAMC ($t_0 =10^{4}$)',
    'samc-1e5-50-slow': 'SAMC ($t_0 =10^{5}$)',
    'samc-1e6-50-slow': 'SAMC ($t_0 =10^{6}$)',
    'samc-1e7-50-slow': 'SAMC ($t_0 =10^{7}$)',
    '1/sqrt(t)': r'$\frac{1}{\sqrt{t}}$',
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
    if fix_legend(m) in _colors:
        return _colors[fix_legend(m)]
    return 'k'

def style_args(method):
    args = {}
    if method != '1/sqrt(t)':
        args = {'label': method}
    if fix_legend(method) in _colors:
        args['color'] = _colors[fix_legend(method)]
    if fix_legend(method) in _linestyles:
        args['linestyle'] = _linestyles[fix_legend(method)]
    args['zorder'] = -legend_order(method)
    if method == '1/sqrt(t)':
        args['zorder'] = -59
        args['linewidth'] = 0.1
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
