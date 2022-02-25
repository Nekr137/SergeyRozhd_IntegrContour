import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def load_columns(fname, delimiter='\t', decimal=','):
    data = pd.read_table(fname, sep=delimiter,header=None, dtype=float, decimal=decimal)
    return [list(col) for col in data.values.T]

def interp(ax, cnt, xdata, ydata):
    xx = np.linspace(xdata[0], xdata[-1], cnt)
    yy = np.interp(xx, xdata, ydata)
    return list(xx),list(yy)

class Perenos:
    def __init__(self, xPos) -> None:
        self._xPos = xPos
        self._depths = []
        self._values = []

    def append(self, depth, value):
        self._depths.append(depth)
        self._values.append(value)

def load_perenos(fname):
    data = load_columns(fname, decimal=',')
    perenos = []
    i = 0
    while i < len(data[0]):
        x = data[0][i]
        perenos.append(Perenos(x))
        while i < len(data[0]) and (data[0][i] - x) < 1e-6:
            depth = data[1][i]
            value = data[2][i]
            perenos[-1].append(depth, value)
            i += 1
    return perenos

def main():
    fig, axs = plt.subplots(nrows=2, ncols=1)
    axs[1].set_xlabel('Distance')
    axs[0].set_ylabel('Depth')
    axs[1].set_ylabel('Depth')

    perenos10 = load_perenos('10.08.copy.perenos.xlsx.dat')
    perenos11 = load_perenos('11.08.copy.perenos.xlsx.dat')
    def show_perenos(ax, perenos):
        for p in perenos:
            for d in p._depths:
                ax.plot([p._xPos], [d], 'k.')
    pnts10 = show_perenos(axs[0], perenos10)
    pnts11 = show_perenos(axs[1], perenos11)

    dno10 = load_columns('dno(10.08.21).bln',delimiter=',', decimal='.')
    dno11 = load_columns('dno(11.08.21).bln',delimiter=',', decimal='.')
    # removing the closing of regions
    dno10 = dno10[0][1:-2], dno10[1][1:-2]
    dno11 = dno11[0][1:-2], dno11[1][1:-2]
    axs[0].plot(dno10[0],dno10[1], 'kx', label='topology')
    axs[1].plot(dno11[0],dno11[1], 'kx', label='topology')

    # interpolating the date
    dno10_interp = interp(axs[0], len(perenos10), dno10[0], dno10[1])
    dno11_interp = interp(axs[1], len(perenos11), dno11[0], dno11[1])
    axs[0].plot(dno10_interp[0], dno10_interp[1], 'r.-', label='interp')
    axs[1].plot(dno11_interp[0], dno11_interp[1], 'r.-', label='interp')

    plt.legend()
    plt.savefig("output_interp")
    pass

if __name__ == '__main__':
    main()