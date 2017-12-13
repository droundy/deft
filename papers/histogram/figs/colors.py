import matplotlib.pyplot as plt

_colors = { 'sad': 'r',
            'sad-tm': 'r',
            'sad3': 'tab:orange',
            'sad3-tm': 'tab:orange',
            'wltmmc-0.8-1e-10': 'k',
            'tmmc': 'b',
            'tmi3': 'g',
            'toe3': 'c',
            'satmmc': 'y',
            'wltmmc-0.8-0.0001': 'tab:purple',
            'wltmmc-1-0.0001': 'tab:pink',
}

_linestyles = {
    'sad-tm': '--',
    'sad3-tm': '--',
}

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
