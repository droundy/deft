import matplotlib.pyplot as plt

_colors = { 'sad': 'r',
            'sad-tm': 'tab:orange',
            'wltmmc-0.8-1e-10': 'k',
            'tmmc': 'b',
            'tmi3': 'g',
            'toe3': 'c',
            'satmmc': 'y',
            'wltmmc-0.8-0.0001': 'tab:purple',
            'wltmmc-1-0.0001': 'tab:pink',
}

def plot(x, y, method=None):
    if method in _colors:
        return plt.plot(x,y, color=_colors[method], label=method)
    return plt.plot(x,y, label=method)

def loglog(x, y, method=None):
    if method in _colors:
        return plt.loglog(x,y, _colors[method], label=method)
    return plt.loglog(x,y, label=method)
