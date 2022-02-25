from pandas import read_table
from numpy import linspace, interp, meshgrid, asfarray, abs, mean
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate

def find_nearest(array, value):
    array = asfarray(array)
    return (abs(array - value)).argmin()

def load_columns(fname, delimiter='\t', decimal=','):
    data = read_table(fname, sep=delimiter,header=None, dtype=float, decimal=decimal)
    return [list(col) for col in data.values.T]

def interp(ax, cnt, xdata, ydata):
    xx = linspace(xdata[0], xdata[-1], cnt)
    yy = interp(xx, xdata, ydata)
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

def show_perenos_pnts(ax, perenos):
    for p in perenos:
        for d in p._depths:
            ax.plot([p._xPos], [d], 'k.')




def build_perenos_data(aPerenos):
    xmin = aPerenos[ 0]._xPos
    xmax = aPerenos[-1]._xPos
    dmin = min([min(p._depths) for p in aPerenos])
    dmax = max([max(p._depths) for p in aPerenos])
    Nx = len(aPerenos)
    Ny = len(aPerenos[0]._depths)
    x = linspace(xmin, xmax, Nx)
    y = linspace(dmin, dmax, Ny)
    X,Y = meshgrid(x,y)
    m = mean([mean(p._values) for p in aPerenos]) # mean values
    D = [[m for i in range(len(x))] for j in range(len(y))]
    for p in aPerenos:
        i = find_nearest(x,p._xPos)
        for idx in range(len(p._depths)):
            value = p._values[idx]
            depth = p._depths[idx]
            j = find_nearest(y, depth)
            D[j][i] = value
    return X,Y,D

def main():
    fig, axs = plt.subplots(nrows=2, ncols=1)
    axs[1].set_xlabel('Distance')
    axs[0].set_ylabel('Depth')
    axs[1].set_ylabel('Depth')

    perenos10 = load_perenos('10.08.copy.perenos.xlsx.dat')
    perenos11 = load_perenos('11.08.copy.perenos.xlsx.dat')
    show_perenos_pnts(axs[0], perenos10)
    show_perenos_pnts(axs[1], perenos11)

    X10,Y10,D10 = build_perenos_data(perenos10)
    X11,Y11,D11 = build_perenos_data(perenos11)

    extent10 = (X10[0][0], X10[0][-1], Y10[0][0], Y10[-1][0])
    extent11 = (X11[0][0], X11[0][-1], Y11[0][0], Y11[-1][0])

    im10 = axs[0].imshow(D10, interpolation='bilinear', origin='lower', cmap=cm.jet, extent=extent10)
    im11 = axs[1].imshow(D11, interpolation='bilinear', origin='lower', cmap=cm.jet, extent=extent11)

    cs10 = axs[0].contour(X10, Y10, D10)
    cs11 = axs[1].contour(X11, Y11, D11)

    axs[0].clabel(cs10, inline=True, fontsize=10)
    axs[1].clabel(cs11, inline=True, fontsize=10)

    plt.savefig("output_interp")
    pass

if __name__ == '__main__':
    main()
