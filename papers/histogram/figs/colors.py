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
            'SAMC ($t_0 =10^{9}$)': 'xkcd:maroon',
            '1/sqrt(t)': 'xkcd:dark gray',
            'converged result': 'xkcd:black',
            r'$\frac{1}{\sqrt{t}}$': '#eeeeee',
}

_linestyles = {
    'SAMC ($t_0 =10^{3}$)': ':',
    'SAMC ($t_0 =10^{4}$)': ':',
    'SAMC ($t_0 =10^{5}$)': ':',
    'SAMC ($t_0 =10^{6}$)': ':',
    'SAMC ($t_0 =10^{7}$)': ':',
    'SAMC ($t_0 =10^{8}$)': ':',
    'SAMC ($t_0 =10^{9}$)': ':',
    'WL': '--',
    '$1/t$-WL': '-.',
    r'$\frac{1}{\sqrt{t}}$': '-',
    'converged result': ':',
}

_legend_order = [
    'bench',

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

    # This is for the Ising paper
    'ising-sad-128',
    'ising-sad-32',
    'ising-wl-128',
    'ising-wl-32',
    'ising-wl-inv-t-128',
    'ising-wl-inv-t-32',
    'ising-samc-1e3-32',
    'ising-samc-1e4-32',
    'ising-samc-1e5-32',
    'ising-samc-1e6-32',
    'ising-samc-1e7-32',
    'ising-samc-1e5-128',
    'ising-samc-1e6-128',
    'ising-samc-1e7-128',
    'ising-samc-1e8-128',
    'ising-samc-1e9-128',

    'lj-sad-31-bin001',
    'lj-sad-31-bin002',

    'lj-wl-31-bin001',
    'lj-wl-31-bin002',
    'lj-wl-31-bin0005',

    'lj-inv-t-wl-31-bin001',
    'lj-inv-t-wl-31-bin002',
    'lj-inv-t-wl-31-bin0005',

    'lj-samc-31-1e5-bin001',
    'lj-samc-31-1e5-bin002',
    'lj-samc-31-1e5-bin0005',

    'lj-samc-31-1e6-bin001',
    'lj-samc-31-1e6-bin002',
    'lj-samc-31-1e6-bin0005',

    'lj-samc-31-1e7-bin001',
    'lj-samc-31-1e7-bin002',
    'lj-samc-31-1e7-bin0005',
]

_legend_label = {
    'sad': 'SAD',
    'sad-256': 'SAD',
    'sad-50': 'SAD',
    'sad-50-slow': 'SAD',
    'wl': 'WL',
    'wl-256': 'WL',
    'wl-50': 'WL',
    'wl-50-slow': 'WL',
    'inv-t-wl': '$1/t$-WL',
    'wl-inv-t-256': '$1/t$-WL',
    'wl-inv-t-50': '$1/t$-WL',
    'wl-inv-t-50-slow': '$1/t$-WL',
    'samc-1e5': 'SAMC ($t_0 =10^{5}$)',
    'samc-1e6': 'SAMC ($t_0 =10^{6}$)',
    'samc-1e7': 'SAMC ($t_0 =10^{7}$)',

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

    # This is for the Ising paper
    'ising-sad-128': 'SAD',
    'ising-sad-32': 'SAD',
    'ising-wl-128': 'WL',
    'ising-wl-32': 'WL',
    'ising-wl-inv-t-128': '$1/t$-WL',
    'ising-wl-inv-t-32': '$1/t$-WL',
    'ising-samc-1e3-32': 'SAMC ($t_0 =10^{3}$)',
    'ising-samc-1e4-32': 'SAMC ($t_0 =10^{4}$)',
    'ising-samc-1e5-32': 'SAMC ($t_0 =10^{5}$)',
    'ising-samc-1e6-32': 'SAMC ($t_0 =10^{6}$)',
    'ising-samc-1e7-32': 'SAMC ($t_0 =10^{7}$)',
    'ising-samc-1e5-128': 'SAMC ($t_0 =10^{5}$)',
    'ising-samc-1e6-128': 'SAMC ($t_0 =10^{6}$)',
    'ising-samc-1e7-128': 'SAMC ($t_0 =10^{7}$)',
    'ising-samc-1e8-128': 'SAMC ($t_0 =10^{8}$)',
    'ising-samc-1e9-128': 'SAMC ($t_0 =10^{9}$)',

    'lj-sad-31-bin001': 'SAD',
    'lj-sad-31-bin002': 'SAD big bin',

    'lj-wl-31-bin001': 'WL',
    'lj-wl-31-bin002': 'WL big bin',
    'lj-wl-31-bin0005': 'WL hires',

    'lj-inv-t-wl-31-bin001': '$1/t$-WL',
    'lj-inv-t-wl-31-bin002': '$1/t$-WL big bin',
    'lj-inv-t-wl-31-bin0005': '$1/t$-WL hires',

    'lj-samc-31-1e5-bin001': 'SAMC ($t_0 =10^{5}$)',
    'lj-samc-31-1e5-bin002': 'SAMC ($t_0 =10^{5}$) big bin',
    'lj-samc-31-1e5-bin0005': 'SAMC ($t_0 =10^{5}$) hires',

    'lj-samc-31-1e6-bin001': 'SAMC ($t_0 =10^{6}$)',
    'lj-samc-31-1e6-bin002': 'SAMC ($t_0 =10^{6}$) big bin',
    'lj-samc-31-1e6-bin0005': 'SAMC ($t_0 =10^{6}$) hires',

    'lj-samc-31-1e7-bin001': 'SAMC ($t_0 =10^{7}$)',
    'lj-samc-31-1e7-bin002': 'SAMC ($t_0 =10^{7}$) big bin',
    'lj-samc-31-1e7-bin0005': 'SAMC ($t_0 =10^{7}$) hires',

    'bench': 'converged result',
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
    args['zorder'] = legend_order(method)
    if fix_legend(method) == 'SAD':
        args['linewidth'] = 3
    if 'SAMC' in fix_legend(method):
        args['linewidth'] = 0.5
    if method == '1/sqrt(t)':
        args['zorder'] = -59
        args['linewidth'] = 1.1
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
